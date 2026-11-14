# HER2 / trastuzumab data exploration

Oslo AI × Antibodies hackathon, 14 November 2026.
Source: `~/hackathon-data/her2-trastuzumab/` (from `s3://apheris-oslo-hackathon/data/`).
Reproduce with `python3 scripts/explore.py` and `python3 scripts/pwm_test.py` (stdlib only).

Both tables come from [Chinery et al.](https://www.biorxiv.org/content/10.1101/2024.03.26.586756v2),
who mutated trastuzumab's CDRH3 (wild type `WGGDGFYAMD`) and measured HER2 binding. Cite them.

---

## 1. The two tables

### SPR affinities — `trastuzumab_her2_spr_all_768.csv`, 768 rows

| `binding_result` | n |
| --- | --- |
| `no_reportable_kd` | 460 |
| `quantified_binding` | 190 |
| `inconclusive` | 104 |
| `not_reported` | 14 |

- CDRH3 lengths: 740 × 10-mer, 20 × 9-mer, 8 × 16-mer. The 16-mers are a different antibody
  (`ARDRIVDDYYYYGMDV`, different VH/VL) — drop them, as the README says.
- 754 unique CDRH3s; 760 of 768 sit on the trastuzumab framework.
- KD range 0.023 – 13,863 nM, median 61, quartiles 7.6 / 193.5. Only 25 variants are sub-nanomolar.

### FACS binding classes — `tz_her2_facs_3classes_524346.csv`, 524,346 rows

- All CDRH3s unique, all exactly 10 residues, single fixed VL.
- Classes: `medium` 190,011 / `high` 172,149 / `low` 162,186.
- Edit distance from wild type is concentrated at 7–8 (311k rows). Only 50 variants at distance 1.

---

## 2. The baseline you have to beat

Trastuzumab appears 8× in the SPR table at **KD 0.029 – 0.178 nM**. Adaptyv will measure it
fresh as the comparison point, so treat ~0.1 nM as the affinity bar.

Only **7 non-reference variants in the entire SPR table** come in below the loosest reference value:

| CDRH3 | KD (nM) | edits from WT |
| --- | --- | --- |
| WGSQGFYAYD | 0.023 | 3 |
| WGSDGWYAMD | 0.061 | 2 |
| WGGQGFYAMD | 0.090 | 1 |
| WGGSGPYLFD | 0.095 | 4 |
| WGGEGFYAMD | 0.117 | 1 |
| WGGDGFYLMD | 0.128 | 1 |
| WGGDGFYEMD | 0.142 | 1 |

---

## 3. The central tension: novelty versus affinity

Affinity decays steeply and monotonically with CDRH3 edit distance:

| edits from WT | n with KD | median KD (nM) | best KD (nM) |
| --- | --- | --- | --- |
| 0 (reference) | 8 | 0.135 | 0.029 |
| 1 | 15 | 4.5 | 0.090 |
| 2 | 14 | 20.5 | 0.061 |
| 3 | 19 | 15.0 | 0.023 |
| 4 | 34 | 60.7 | 0.095 |
| 5 | 32 | 92.7 | 0.471 |
| 6 | 29 | 93.0 | 0.185 |
| 7 | 18 | 149.6 | 2.9 |
| 8 | 10 | 248.5 | 30.0 |
| 9–10 | 8 | ~200 | 66 |

The optimisation-track novelty filter requires **under 70% CDRH3 identity to trastuzumab**, i.e.
**≥4 mutations on a 10-mer**. At ≥4 edits the measured median KD is ~116 nM — roughly 1000× worse
than trastuzumab.

**The window is thin but not empty.** Of the 131 novelty-passing SPR variants with a KD,
**9 are under 1 nM (6.9%)** and the best, `WGGSGPYLFD`, reaches 0.095 nM at 4 edits — competitive
with the reference. That 6.9% is the number any ranking method has to beat: pick at random inside
the novelty window and roughly 1 in 14 is sub-nanomolar.

---

## 4. Per-position preferences (FACS, high vs low)

Log-ratio enrichment per position, `high` against `low`:

| pos | WT | favoured in `high` | disfavoured |
| --- | --- | --- | --- |
| 1 | W | N 2.3× Y 1.8× H 1.6× | I 0.06× A 0.08× |
| 2 | G | A 2.5× G 2.3× D 2.1× | K 0.23× W 0.23× |
| 3 | G | **A 10.7×** G 3.4× R 2.7× | H 0.15× C 0.31× |
| 4 | D | G 4.5× P 2.2× A 2.2× | H 0.18× E 0.26× |
| 5 | G | **R 13.6×** P 2.7× M 2.3× | K 0.07× Y 0.14× |
| 6 | F | G 3.4× D 2.3× P 2.3× | K 0.26× Q 0.27× |
| 7 | Y | W 1.9× G 1.8× D 1.6× | T 0.29× K 0.35× |
| 8 | A | T 4.3× S 3.3× K 3.2× | L 0.27× F 0.29× |
| 9 | M | A 1.6× V 1.5× L 1.3× | W 0.39× H 0.44× |
| 10 | D | D 5.9× N 3.0× Y 1.9× | W 0.20× C 0.20× |

Positions 3, 5, 8 and 10 carry real signal; 1 and 9 are near-indifferent. Lysine is disfavoured
almost everywhere. Useful as a **generation prior** — see the caveat below before using it to rank.

---

## 5. Two negative results that should change your pipeline

### 5.1 FACS class does not predict affinity

37 CDRH3s appear in both tables; 27 have an SPR outcome.

| FACS class | n with KD | median KD (nM) | range (nM) |
| --- | --- | --- | --- |
| high | 16 | 8.0 | 0.09 – 174 |
| medium | 6 | 4.6 | 0.47 – 139 |
| low | 5 | 15.4 | 0.44 – 94 |

`medium` has a *better* median than `high`, the ranges overlap almost entirely, and 4 of 21
FACS-`high` variants gave no reportable KD at all. The README states the classes report whether a
variant binds, not how tightly — this is the size of that gap. (n is 5–16 per class, so this bounds
the usefulness of FACS class as an affinity label rather than measuring the gap precisely.)

### 5.2 A FACS-derived PWM does not rank affinity either

I built a position-weight matrix of log-odds (`high` vs `low`) from all 524k FACS rows and scored the
179 trastuzumab-framework 10-mer SPR variants that have a KD (reference rows excluded).

| predictor | Spearman ρ vs log₁₀ KD |
| --- | --- |
| PWM log-odds | **−0.116** |
| edit distance from WT | +0.412 |

The PWM is essentially uninformative, and what little it has is carried by edit distance — within
edit-distance strata ρ is near zero or the wrong sign (+0.22, +0.10, +0.23, +0.12, +0.17, −0.27,
+0.23 for distances 1–7). Its top decile has a *worse* median KD (60.9 nM) than you would get at
random, and it ranks the 7 best-known binders at positions 48, 53, 54, 82, 92, 102 and 111 of 179 —
the middle of the pack. Within the novelty window (≥4 edits) ρ = −0.043.

**Conclusion: do not rank designs by a FACS-trained classifier score.** It is a reasonable prior for
*generating* plausible loops, and a reasonable binder/non-binder gate, but it carries no affinity
ranking signal. Edit distance correlates with affinity but is unusable as a ranker, because the
novelty filter forces every submission to ≥4 edits.

---

## 6. What the data supports for ranking designs

Ranked by how much evidence there is behind each:

1. **Rank on the 190 SPR KDs, not on FACS.** This is the only affinity ground truth in the
   dataset. It is small, so use heavy regularisation and nested cross-validation, and treat the
   output as a coarse filter (3–5 bins) rather than a fine ordering. Hold out the 9 sub-nanomolar
   novelty-passing variants as your acceptance test: a ranker that cannot surface those is not
   ranking affinity.
2. **Use the 460 `no_reportable_kd` rows as negatives.** Expressed, measured, did not bind — the
   most valuable label in the table, and far more numerous than the positives.
3. **Calibrate structural confidence before trusting it.** `~/hackathon-data/sigma2-ab-ag-champloo/`
   and `~/hackathon-data/zenodo-18390239-ab-ag-complexes/` hold 454 SPR-measured non-binders (21
   scoring AF3 ipTM ≥ 0.7) and 1.8M predictions over 106 cognate / 11,342 non-cognate pairs. Fit
   your ipTM / min_ipae / ipSAE threshold on that, rather than assuming one. Rank within a single
   folding model only — scores are not comparable across models.
4. **Gate on `abguard` first, not last.** It decides your track and never reclassifies. Most
   conservative and balanced `her2-2026-11-07` output fails novelty (6% and 0% pass rates per the
   task PDF); `her2-novel-2026-11-09` passes 86%.
5. **Screen for expression before submitting.** The task PDF reports that Adaptyv expression
   success improved once teams scored candidates with Aggrescan 3D, SoluProt, NetSolP and
   soluble-weights ProteinMPNN. A design that does not express scores zero on all three judged axes.
6. **Do not neglect thermostability.** It is one of the three judged properties and nothing in these
   two tables measures it. Without a proxy for it you are capped at A-tier.

### Suggested ranking function

```
eligible  = abguard pass (hard gate, first)
score     = w1 · SPR-trained affinity model (binned)
          + w2 · calibrated interface confidence (min_ipae / ipSAE, one folding model)
          + w3 · expression/solubility panel
          + w4 · thermostability proxy
```

Keep the weights explicit and write down why you chose them — judgment is 40% of the pitch score,
and "we ranked by ipTM" is the answer the judges will have heard most often.

---

## 7. Caveats

- Figures and numbers from the task PDF were extracted with a custom stdlib parser (no poppler on
  the VM); body text was clean, figure wording is approximate.
- The FACS↔SPR comparison rests on 27 overlapping measurements. It is enough to show FACS class is
  not an affinity label, not enough to quantify the relationship.
- The PWM is a linear, position-independent model trained on class labels. A richer model on the
  FACS data might extract affinity signal the PWM cannot — but it would still be learning from
  binding classes that themselves do not track KD, so I would not expect much.
- `edit_distance` in the FACS table is Hamming distance on a fixed 10-mer; CDRH3 novelty under
  `abguard` is sequence identity, which is the same thing here but will not be for 9-mers or
  indel-containing designs.
