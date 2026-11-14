#!/usr/bin/env python3
"""
Train a FACS binding-class gate for trastuzumab CDRH3 variants.

Task: binary, `high` vs {`mid`, `low`} -> one probability P(high).

Trains four models plus a no-model floor, all on the same clonotype split, and
reports every one on val AND test. The winner is chosen on **val** (average
precision); test is read once and reported. Choosing on test would make the
headline number optimistic.

What this is NOT: an affinity predictor. The training signal is a FACS sort gate,
which scores Spearman -0.116 against log KD on this project's own data (see
../her2-data-exploration.md). P(high) is a predicted binding class, not a KD.

Usage:
    python train.py --quick     # debug on train1120_val240, seconds
    python train.py             # full 0.7_0.15_0.15 split
"""
import argparse
import copy
import json
import os
import time

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score
from sklearn.neural_network import MLPClassifier

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}
TZ_CDRH3 = "WGGDGFYAMD"      # trastuzumab wild type
CDRH3_LEN = 10
SEED = 0

SPLIT_BASE = (
    "/home/hackathon/hackathon-data/zenodo-10831512-trastuzumab-ngs/ngs_unz/data/"
    "affinity_data/her2/her2_aff_large/processed/remove_overlap/clonotype_split"
)
HERE = os.path.dirname(os.path.abspath(__file__))


# --------------------------------------------------------------------------- #
# features
# --------------------------------------------------------------------------- #
def onehot(seqs):
    """(n, 200) one-hot of a fixed-length-10 CDRH3."""
    X = np.zeros((len(seqs), CDRH3_LEN, 20), dtype=np.float32)
    for r, s in enumerate(seqs):
        for p, a in enumerate(s):
            j = AA_IDX.get(a)
            if j is not None:
                X[r, p, j] = 1.0
    return X.reshape(len(seqs), CDRH3_LEN * 20)


def edit_distance(seqs):
    """Mutations vs trastuzumab's CDRH3, scaled to [0,1] so it does not dominate
    the regularisation of the linear/MLP models (one-hot columns are 0/1)."""
    d = [sum(1 for a, b in zip(s, TZ_CDRH3) if a != b) for s in seqs]
    return np.asarray(d, dtype=np.float32).reshape(-1, 1) / CDRH3_LEN


def featurise(seqs):
    """201 features: 200 one-hot + 1 scaled edit distance."""
    return np.hstack([onehot(seqs), edit_distance(seqs)])


FEATURE_SPEC = {
    "cdrh3_len": CDRH3_LEN,
    "alphabet": AA,
    "layout": "onehot(10x20) flattened, then edit_distance_to_trastuzumab/10",
    "n_features": CDRH3_LEN * 20 + 1,
    "reference_cdrh3": TZ_CDRH3,
}


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #
def load_splits(split_dir):
    out = {}
    for name in ("train", "val", "test"):
        df = pd.read_csv(os.path.join(split_dir, f"{name}.csv"))
        before = len(df)
        df = df[df["seq"].str.len() == CDRH3_LEN].copy()
        if len(df) != before:
            print(f"  {name}: dropped {before - len(df)} rows of length != {CDRH3_LEN}")
        out[name] = df
    return out


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #
def score(y, p):
    return {"roc_auc": float(roc_auc_score(y, p)),
            "avg_precision": float(average_precision_score(y, p))}


def best_f1_threshold(y, p):
    """Threshold maximising F1 on the given (validation) split."""
    prec, rec, thr = precision_recall_curve(y, p)
    # precision_recall_curve returns len(thr) == len(prec) - 1
    prec, rec = prec[:-1], rec[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        f1 = np.where((prec + rec) > 0, 2 * prec * rec / (prec + rec), 0.0)
    i = int(np.nanargmax(f1))
    return float(thr[i]), float(prec[i]), float(rec[i]), float(f1[i])


def operating_points(y, p, thresholds=(0.3, 0.4, 0.5, 0.6, 0.7)):
    """Precision / recall / pass-rate at fixed thresholds, so a permissive gate
    can be chosen deliberately rather than by accepting 0.5."""
    rows = []
    for t in thresholds:
        pred = p >= t
        tp = int(np.sum(pred & (y == 1)))
        fp = int(np.sum(pred & (y == 0)))
        fn = int(np.sum(~pred & (y == 1)))
        prec = tp / (tp + fp) if tp + fp else float("nan")
        rec = tp / (tp + fn) if tp + fn else float("nan")
        rows.append({"threshold": t, "precision": prec, "recall": rec,
                     "pass_rate": float(np.mean(pred))})
    return rows


# --------------------------------------------------------------------------- #
# models
# --------------------------------------------------------------------------- #
def fit_mlp(Xtr, ytr, Xva, yva, max_epochs=30, patience=3):
    """
    MLP with early stopping on the REAL validation split.

    sklearn's MLPClassifier(early_stopping=True) validates on an internal random
    split of the training data. Under a clonotype split that leaks one-mutation
    neighbours across the stopping criterion, so we drive the epochs ourselves
    with partial_fit (one epoch per call) and stop on val average precision.
    """
    clf = MLPClassifier(hidden_layer_sizes=(256, 128), alpha=1e-4,
                        batch_size=512, learning_rate_init=1e-3,
                        random_state=SEED)
    rng = np.random.default_rng(SEED)
    best_ap, best_model, best_epoch, bad = -np.inf, None, 0, 0
    classes = np.array([0, 1])
    for epoch in range(1, max_epochs + 1):
        idx = rng.permutation(len(Xtr))
        clf.partial_fit(Xtr[idx], ytr[idx], classes=classes)
        ap = average_precision_score(yva, clf.predict_proba(Xva)[:, 1])
        if ap > best_ap + 1e-5:
            best_ap, best_model, best_epoch, bad = ap, copy.deepcopy(clf), epoch, 0
        else:
            bad += 1
            if bad >= patience:
                break
    print(f"      MLP: best val AP {best_ap:.4f} at epoch {best_epoch} "
          f"(ran {epoch})")
    return best_model, best_epoch


def fit_lightgbm(Xtr, ytr, Xva, yva):
    import lightgbm as lgb
    clf = lgb.LGBMClassifier(n_estimators=300, learning_rate=0.1, num_leaves=63,
                             objective="binary", n_jobs=-1, verbose=-1,
                             random_state=SEED)
    clf.fit(Xtr, ytr, eval_X=Xva, eval_y=yva, eval_metric="average_precision",
            callbacks=[lgb.early_stopping(30, verbose=False)])
    print(f"      LightGBM: best iteration {clf.best_iteration_} of 300")
    return clf


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="debug on train1120_val240 instead of the full split")
    ap.add_argument("--split-dir", default=None, help="override the split directory")
    ap.add_argument("--out-dir", default=os.path.join(HERE, "model"))
    args = ap.parse_args()

    split_dir = args.split_dir or os.path.join(
        SPLIT_BASE, "train1120_val240" if args.quick else "0.7_0.15_0.15")
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"split: {split_dir}")
    d = load_splits(split_dir)
    for k, v in d.items():
        pos = int(v["label"].sum())
        print(f"  {k:5} n={len(v):>7}  positives={pos:>7} ({100 * pos / len(v):.2f}%)  "
              f"clusters={v['cluster_centre'].nunique()}")

    Xtr, ytr = featurise(d["train"]["seq"].tolist()), d["train"]["label"].to_numpy()
    Xva, yva = featurise(d["val"]["seq"].tolist()), d["val"]["label"].to_numpy()
    Xte, yte = featurise(d["test"]["seq"].tolist()), d["test"]["label"].to_numpy()

    # sanity: one-hot part must have exactly one 1 per position
    rs = Xtr[:, :200].sum(axis=1)
    assert np.allclose(rs, CDRH3_LEN), f"one-hot row sums wrong: {rs.min()}..{rs.max()}"
    print(f"  one-hot check ok (row sums == {CDRH3_LEN}); features = {Xtr.shape[1]}")

    results, fitted = [], {}

    def record(name, p_va, p_te, secs, model=None, note=""):
        row = {"model": name, "val": score(yva, p_va), "test": score(yte, p_te),
               "train_seconds": round(secs, 1), "note": note}
        results.append(row)
        fitted[name] = model
        print(f"  {name:26} val AP={row['val']['avg_precision']:.4f} "
              f"AUC={row['val']['roc_auc']:.4f} | "
              f"test AP={row['test']['avg_precision']:.4f} "
              f"AUC={row['test']['roc_auc']:.4f}  [{secs:.1f}s]")

    print("\nmodels (val chooses the winner; test is reported once):")

    # 0 - floor: fewer mutations from trastuzumab => more likely `high`
    record("0 edit distance (floor)",
           -Xva[:, -1], -Xte[:, -1], 0.0, None,
           "no model; score = -edit_distance")

    # 1 - logistic on edit distance only
    t0 = time.time()
    m1 = LogisticRegression(max_iter=1000).fit(Xtr[:, -1:], ytr)
    record("1 logistic (edit dist)",
           m1.predict_proba(Xva[:, -1:])[:, 1],
           m1.predict_proba(Xte[:, -1:])[:, 1], time.time() - t0, m1)

    # 2 - logistic on one-hot + edit distance
    t0 = time.time()
    m2 = LogisticRegression(max_iter=400).fit(Xtr, ytr)
    record("2 logistic (one-hot)",
           m2.predict_proba(Xva)[:, 1], m2.predict_proba(Xte)[:, 1],
           time.time() - t0, m2)

    # 3 - MLP, early stopped on the real val split
    t0 = time.time()
    m3, ep = fit_mlp(Xtr, ytr, Xva, yva)
    record("3 MLP 256-128",
           m3.predict_proba(Xva)[:, 1], m3.predict_proba(Xte)[:, 1],
           time.time() - t0, m3, f"early stopped at epoch {ep} on val AP")

    # 4 - LightGBM
    try:
        t0 = time.time()
        m4 = fit_lightgbm(Xtr, ytr, Xva, yva)
        record("4 LightGBM",
               m4.predict_proba(Xva)[:, 1], m4.predict_proba(Xte)[:, 1],
               time.time() - t0, m4, f"best_iteration={m4.best_iteration_}")
    except ImportError:
        print("  LightGBM unavailable, skipped")

    # ---- pick the winner on VAL -------------------------------------------- #
    trainable = [r for r in results if fitted.get(r["model"]) is not None]
    winner = max(trainable, key=lambda r: r["val"]["avg_precision"])
    floor = results[0]
    by_test = max(trainable, key=lambda r: r["test"]["avg_precision"])

    print("\n" + "=" * 94)
    print(f"{'model':26} {'val AP':>9} {'val AUC':>9} {'test AP':>9} {'test AUC':>9} {'secs':>7}")
    print("-" * 94)
    for r in results:
        print(f"{r['model']:26} {r['val']['avg_precision']:>9.4f} "
              f"{r['val']['roc_auc']:>9.4f} {r['test']['avg_precision']:>9.4f} "
              f"{r['test']['roc_auc']:>9.4f} {r['train_seconds']:>7}")
    print("=" * 94)

    print(f"\nwinner on val: {winner['model']}")
    if by_test["model"] != winner["model"]:
        print(f"NOTE: test would have preferred '{by_test['model']}'. Val and test "
              f"disagree -- reporting the val-selected model, as planned. Val has "
              f"fewer clonotype clusters, so its ranking is the noisier one.")
    if winner["val"]["avg_precision"] <= floor["val"]["avg_precision"]:
        print("WARNING: no model beats the edit-distance floor on val. "
              "Prefer the floor and do not ship a model.")

    # ---- threshold from VAL ------------------------------------------------ #
    p_va = (fitted[winner["model"]].predict_proba(Xva)[:, 1])
    thr, prec, rec, f1 = best_f1_threshold(yva, p_va)
    print(f"\nthreshold from val (max F1): {thr:.4f}  "
          f"precision={prec:.3f} recall={rec:.3f} F1={f1:.3f}")
    ops_va = operating_points(yva, p_va)
    print("\nval operating points (pick a permissive gate deliberately):")
    print(f"  {'thr':>5} {'precision':>10} {'recall':>8} {'pass rate':>10}")
    for r in ops_va:
        print(f"  {r['threshold']:>5.2f} {r['precision']:>10.3f} "
              f"{r['recall']:>8.3f} {r['pass_rate']:>10.3f}")

    p_te = fitted[winner["model"]].predict_proba(Xte)[:, 1]
    ops_te = operating_points(yte, p_te, thresholds=(thr,))
    print(f"\nat that threshold on test: precision={ops_te[0]['precision']:.3f} "
          f"recall={ops_te[0]['recall']:.3f} pass_rate={ops_te[0]['pass_rate']:.3f}")

    # ---- persist ----------------------------------------------------------- #
    meta = {
        "task": "binary high vs {mid, low} FACS binding class",
        "warning": "P(high) is a predicted FACS binding class, NOT an affinity/KD estimate",
        "split_dir": split_dir,
        "quick": bool(args.quick),
        "seed": SEED,
        "feature_spec": FEATURE_SPEC,
        "n_train": int(len(ytr)), "n_val": int(len(yva)), "n_test": int(len(yte)),
        "positive_rate_test": float(yte.mean()),
        "selected_on": "val average_precision",
        "winner": winner["model"],
        "threshold": thr,
        "threshold_source": "max F1 on val",
        "val_operating_points": ops_va,
        "test_at_threshold": ops_te[0],
        "results": results,
    }
    with open(os.path.join(args.out_dir, "metrics.json"), "w") as fh:
        json.dump(meta, fh, indent=2)

    joblib.dump({"model": fitted[winner["model"]], "model_name": winner["model"],
                 "threshold": thr, "feature_spec": FEATURE_SPEC, "meta": meta},
                os.path.join(args.out_dir, "gate.joblib"))

    print(f"\nwrote {os.path.join(args.out_dir, 'gate.joblib')}")
    print(f"wrote {os.path.join(args.out_dir, 'metrics.json')}")


if __name__ == "__main__":
    main()
