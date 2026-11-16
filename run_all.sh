#!/usr/bin/env bash
# Reproduces every result in the paper from the two public sources. Run from the repository root.
# R_NULL (default 500) and B_DMACS (default 300) set the simulation sizes used in the paper.
set -euo pipefail
R_NULL=${R_NULL:-500}; B_DMACS=${B_DMACS:-300}; SEED=20261113
mkdir -p derived out fig
bash scripts/00_download_data.sh
HF=data/hf/contents.parquet; EP=data/epoch/extracted
python3 scripts/01_prep_hf.py  $HF derived                        # hf_models_v1.csv, hf_pairs_v1.csv
python3 scripts/02_prep_logit.py $HF derived/hf_models_v1.csv derived   # hf_models_v2.csv (logit scale)
python3 scripts/03_reliability.py $HF out/retest_logit_v1.csv out/precision_diffs_v1.csv
Rscript scripts/04_descriptives.R
Rscript scripts/05_select_raw.R   derived/hf_models_v1.csv out     # raw-scale selection (sensitivity)
Rscript scripts/06_select_logit.R derived/hf_models_v2.csv out     # M4 selection and holdout confirmation
M4=out/m4_syntax_logit_v1.txt
Rscript scripts/07_invariance_role.R   derived/hf_models_v2.csv out $M4
Rscript scripts/08_invariance_cohort.R derived/hf_models_v2.csv out $M4
Rscript scripts/09_mcv_partial.R derived/hf_models_v2.csv out/retest_logit_v1.csv out $M4
Rscript scripts/12_kenny_dimensionality.R derived/hf_models_v2.csv out/g_kenny_v1.txt   # writes out/m5_syntax_logit_v1.txt
Rscript scripts/13_twofactor.R                                     # writes out/m2f_syntax_v1.txt
D=derived/hf_models_v2.csv; N=scripts/10_invariance_null.R
Rscript $N $D $M4 role   $R_NULL out/e_mcnull_role_exact_v1.txt   $SEED exact
Rscript $N $D $M4 cohort $R_NULL out/e_mcnull_cohort_exact_v1.txt $SEED exact
Rscript $N $D $M4 role   $R_NULL out/e_mcnull_role_ext_v1.txt     $SEED ext
Rscript $N $D $M4 cohort $R_NULL out/e_mcnull_cohort_ext_v1.txt   $SEED ext
Rscript $N $D $M4 role   $R_NULL out/e_mcnull_role_misfit_v1.txt  $SEED misfit
Rscript $N $D $M4 cohort $R_NULL out/e_mcnull_cohort_misfit_v1.txt $SEED misfit   # construction check; reported as failed
Rscript $N $D out/m2f_syntax_v1.txt role $(( R_NULL * 3 / 5 )) out/e_mcnull_role_2f_misfit_v1.txt $SEED misfit
Rscript scripts/11_dmacs.R $D out/f_dmacs_v1.txt $B_DMACS
Rscript scripts/14_mimic_size_sensitivity.R $D out/h_ritchie_v1.txt
Rscript scripts/15_pairs_selection_ptonly.R
Rscript scripts/16_mcv_strict.R > out/m_mcv_strict_v1.txt
Rscript scripts/17_robust_exact.R > out/o_robust_exact_v1.txt
python3 scripts/18_epoch_long.py $EP derived
python3 scripts/19_epoch_prep.py $EP derived
python3 scripts/19_epoch_prep.py $EP derived median                  # median-collapse robustness file
Rscript scripts/20_epoch.R
Rscript scripts/21_epoch_noclip.R
Rscript scripts/22_g_change_controls.R
Rscript scripts/23_q_calibration.R > out/p_qcal_v1.txt
Rscript scripts/24_placebo.R
python3 scripts/27_prep_extra.py $HF derived                       # precision, chat template, weight type
Rscript scripts/28_review_tests.R                                  # writes derived/hf_models_v2_notemplate.csv
Rscript $N $D $M4 role $(( R_NULL * 2 / 5 )) out/e_mcnull_role_org_exact_v1.txt $SEED exact org
Rscript $N derived/hf_models_v2_notemplate.csv $M4 role $(( R_NULL * 2 / 5 )) out/e_mcnull_role_notemplate_exact_v1.txt $SEED exact
Rscript scripts/29_epoch_tests.R $(( R_NULL * 2 )) $R_NULL
python3 scripts/30_prep_norm.py $HF derived/hf_models_v1.csv derived  # normalised scores for the clipping-constant test
Rscript scripts/30_review_tests2.R                                 # writes the matched and equal-level subsets
Rscript $N derived/hf_models_v2_matched_role.csv $M4 role $(( R_NULL * 2 / 5 )) out/e_mcnull_role_matched_exact_v1.txt $SEED exact
Rscript $N derived/hf_models_v2_matched_cohort.csv $M4 cohort $(( R_NULL * 2 / 5 )) out/e_mcnull_cohort_matched_exact_v1.txt $SEED exact
Rscript $N derived/hf_models_v2_cohort_2223_25.csv $M4 cohort $(( R_NULL * 2 / 5 )) out/e_mcnull_cohort_2223v25_exact_v1.txt $SEED exact
Rscript scripts/31_epoch_missing.R
Rscript scripts/32_review_tests3.R                                 # selection of descendants; Holm adjustment (reads the null outputs)
Rscript scripts/33_review_tests4.R all out/w_review_tests4_v1.txt   # bootstrap coverage calibration; evaluation period and drift; writes the within-period subsets
Rscript $N derived/hf_models_v2_period_early.csv $M4 role $(( R_NULL * 2 / 5 )) out/e_mcnull_role_early_exact_v1.txt $SEED exact
Rscript $N derived/hf_models_v2_period_late.csv $M4 role $(( R_NULL * 2 / 5 )) out/e_mcnull_role_late_exact_v1.txt $SEED exact
Rscript $N derived/hf_models_v2_cohort_early.csv $M4 cohort $(( R_NULL * 2 / 5 )) out/e_mcnull_cohort_early_exact_v1.txt $SEED exact
Rscript $N derived/hf_models_v2_cohort_late.csv $M4 cohort $(( R_NULL * 2 / 5 )) out/e_mcnull_cohort_late_exact_v1.txt $SEED exact
python3 scripts/25_fig1.py
python3 scripts/26_fig2.py
sha256sum out/* fig/* > out/SHA256SUMS
