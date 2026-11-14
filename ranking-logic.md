# Design ranking logic

`scripts/score_designs_ziyu.py` → `scored_designs.csv`, ranking 142 designs.

## Three inputs

| File | Score | What it measures |
|---|---|---|
| `iptm_binder_antigen.csv` | `iptm_binder_antigen_max` | VH-to-HER2 interface |
| `passing.with_p_high.csv` | `p_high` | Sequence classifier |
| `designs_combined.csv` | liability, ddG, CDRH3 confidence, solubility | How manufacturable the antibody is |

The three are uncorrelated, so they are independent evidence.

## Four steps

**1. Developability gates** → 118
Hard exclusion, not scored: `tier`/`tier_apo` = pass, no major flags, `liability_score` ≤ 4, ddG per substitution ≤ 0.60. A molecule that is hard to make is unusable however well it binds.

**2. Interface gate** → 27
Interface score must exceed **0.231**, which is what trastuzumab itself scores under this pipeline. The absolute ipTM scale is compressed here — SPR-validated sub-nanomolar binders score only 0.181–0.282 — so the pass mark is 0.231, **not** the 0.8 usually quoted.

**3. Remove redundancy** → 20
Keep the best design per sequence; group by first residue + last four, cap each CDRH3 family at 2.

**4. Rank by composite score**
Each component becomes a percentile within the full 142 (raw units are incommensurable), then:

```
score = 0.50 × interface        iptm_binder_antigen_max
      + 0.20 × CDRH3 confidence cdrh3_plddt_mean
      + 0.15 × thermostability  thermo_ddg_per_subst  (lower is better)
      + 0.15 × solubility       mpnn_delta_vs_wt      (lower is better)
```

Tiers A/B/C (> 0.40 / > 0.282 / > 0.231) annotate confidence only; they do not affect the ordering.

## Why p_high is not scored

Spearman ρ = **−0.067** against the interface score — no relationship. Adding it would only inject noise. It is carried as an annotation and a last-resort tiebreak. Only 2 of the 20 selected designs pass its gate; that disagreement is worth reporting, not averaging away.

## Three limitations

- **Epitope unverified**: `target_aligned_binder_rmsd` is empty (target not fixed), so this answers "does the interface look right", not "does it hit the epitope".
- **The 0.231 baseline** comes from `reference_scores.csv`, not from the three inputs. It is the constant `REF_TRASTUZUMAB`, overridable with `--ref-iptm`.
- **142 ≠ 200**: of the 58 missing, 39 have `p_high` ≥ 0.9, so the gap is biased against high-scoring designs.

Every threshold is a CLI flag; `--keep-all` writes all 142 rows with the gate that removed each.
