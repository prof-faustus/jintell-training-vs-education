# Training versus Education in Machines — analysis code

Code that reproduces every number, table and figure in:

Wright, C. S. (2026). *Training versus Education in Machines: Fine-Tuning and Generational Gains in Large Language Models Are Not Measurement-Invariant Gains in g* Manuscript prepared for submission to *Journal of Intelligence* (MDPI).

The analysis uses only two public secondary sources. There are no human participants.

## Data

| Source | What is used | Licence | In this repository |
|---|---|---|---|
| Hugging Face Open LLM Leaderboard `open-llm-leaderboard/contents`, revision `9c09a7cae43334062a82cb164f2ef255013dafa2` (last modified 2025-03-20), file `data/train-00000-of-00001.parquet` | Six benchmark scores (IFEval, BBH, MATH Lvl 5, GPQA, MuSR, MMLU-Pro) with model type, base model, generation, precision and upload date; 4,576 rows, giving 4,496 unique unflagged models | No licence declared on the dataset card | **Not redistributed.** `scripts/00_download_data.sh` downloads the pinned revision and checks SHA-256 `f13a7061e63133bc8745d89f8c7a56148584940ede22d3b27e21cad1e078dca5`. |
| Epoch AI Benchmarking Hub, `benchmark_data.zip` | Frontier-model scores on GPQA Diamond, OTIS Mock AIME 2024–2025, ARC-AGI, ARC-AGI-2 and DTBench, with each benchmark's random baseline and ceiling | CC BY 4.0 | Included as `data/epoch/benchmark_data.zip` (SHA-256 `7fd923b837532422a0a35aabec72c753b6f554fd425221fa2bb3aeb9948ffe34`), unmodified. Attribution: Epoch AI, 'Capabilities & benchmarking', published online at epoch.ai, retrieved from https://epoch.ai/benchmarks. |

If the Hugging Face file at the pinned revision ever changes, the hash check fails and the pipeline stops.

## Requirements

- R 4.3.3 with lavaan 0.6-17, psych 2.4.1 and MASS 7.3-60.0.1
- Python 3.13 with the packages in `requirements.txt`
- bash, curl, unzip and sha256sum

## Running

```bash
bash run_all.sh                         # full run with the simulation sizes used in the paper
R_NULL=20 B_DMACS=20 bash run_all.sh    # quick check with reduced simulation sizes
```

Run it from the repository root. Outputs go to `derived/` (prepared data), `out/` (results) and `fig/` (figures), and `out/SHA256SUMS` lists their hashes. All random steps use seed 20261113. The full run is dominated by the simulated invariance nulls (`10_invariance_null.R`: 500 replicates for each of six M4 nulls, 300 for the two-factor null, 200 for the organisation-clustered and no-template nulls) and the 1000-permutation frontier tests (`29_epoch_tests.R`).

## Scripts and where their results appear

| Script | Output | Paper |
|---|---|---|
| `00_download_data.sh` | `data/hf/contents.parquet`, `data/epoch/extracted/` | Data |
| `01_prep_hf.py` | `derived/hf_models_v1.csv`, `derived/hf_pairs_v1.csv` | Sample construction |
| `02_prep_logit.py` | `derived/hf_models_v2.csv` | Logit scale and clipping |
| `03_reliability.py` | `out/retest_logit_v1.csv`, `out/precision_diffs_v1.csv` | Precision-duplicate reliability; negative control |
| `04_descriptives.R` | `out/descriptives_v1.txt` | Correlations, Mardia, lineage ICCs |
| `05_select_raw.R`, `06_select_logit.R` | `out/b_select*_v1.txt`, `out/m4_syntax*_v1.txt` | Model selection on the discovery half, holdout confirmation |
| `07_invariance_role.R`, `08_invariance_cohort.R` | `out/c_invariance*_v1.txt` | Invariance sequences (role, cohort) |
| `09_mcv_partial.R` | `out/d_mcv_partial_v1.txt` | Method of correlated vectors; partial invariance |
| `10_invariance_null.R` | `out/e_mcnull_*_v1.txt` and `_reps.csv` | Simulated Δfit nulls (exact, ext and misfit constructions) |
| `11_dmacs.R` | `out/f_dmacs_v1.txt` | dMACS effect sizes under six referents |
| `12_kenny_dimensionality.R` | `out/g_kenny_v1.txt`, `out/m5_syntax_logit_v1.txt` | Dimensionality; M5 generating model |
| `13_twofactor.R` | `out/g2_twofactor*_v1.txt`, `out/m2f_syntax_v1.txt` | Two-factor alternative and its cohort inadmissibility |
| `14_mimic_size_sensitivity.R` | `out/h_ritchie_v1.txt` | MIMIC models A/B/C |
| `15_pairs_selection_ptonly.R` | `out/j_extra_v1.txt` | Strict pairs; base selection; pretrained-only role contrast |
| `16_mcv_strict.R` | `out/m_mcv_strict_v1.txt` | MCV on strict pairs |
| `17_robust_exact.R` | `out/o_robust_exact_v1.txt` | Robustness of the exact null |
| `18_epoch_long.py`, `19_epoch_prep.py`, `20_epoch.R`, `21_epoch_noclip.R` | `derived/long_v1.csv`, `derived/epoch_wide_v1.csv`, `out/i_epoch*_v1.txt` | Epoch FIML replication |
| `22_g_change_controls.R` | `out/p_strong_v1.txt` | Within-lineage g change, Q test, TOST, negative and positive controls, depth |
| `23_q_calibration.R` | `out/p_qcal_v1.txt` | Empirical calibration of Q |
| `24_placebo.R` | `out/q_placebo_v1.txt` | Placebo splits |
| `27_prep_extra.py` | `derived/hf_extra_v1.csv` | Evaluation precision, chat-template flag, weight type |
| `28_review_tests.R` | `out/r_review_tests_v1.txt` | Type-label check; evaluation format; matched pairs; organisation-level dependence; linking constraints (Table 5); IFEval contrast; base-selection coverage |
| `10_invariance_null.R` (organisation clusters; no-template sample) | `out/e_mcnull_role_org_exact_v1.txt`, `out/e_mcnull_role_notemplate_exact_v1.txt` | Training-role nulls under organisation clustering and without chat-template evaluations |
| `29_epoch_tests.R` | `out/s_epoch_tests_v1.txt` | Frontier replication: configural fit, permutation tests, median collapse |
| `25_fig1.py`, `26_fig2.py` | `fig/fig1_v1.pdf`, `fig/fig2_v1.pdf` | Figures 1 and 2 |

## Reference results

`reference_results/` holds the output of a full run, with SHA-256 hashes in `reference_results/SHA256SUMS`, so a rerun can be compared file by file. The non-random outputs were regenerated from a clean copy of this repository and are byte-identical across runs. For the simulations, a reduced run (`R_NULL=20`) reproduced exactly the first replicates of the seven nulls that predate the review tests, because each null draws from a seeded random stream. The two later nulls (organisation clusters; no chat template), the review tests and the frontier permutation tests were run at full size with the scripts in this repository. The files contain summary statistics only and no rows of the Hugging Face table.

## Licence

The code is released under the MIT licence (`LICENSE`). The licence does not cover the data. The Epoch AI file stays under CC BY 4.0. The Hugging Face data is downloaded from its source and is not part of this repository.
