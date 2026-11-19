# Supporting Information

## MLPATH: Energy-and-force transition-state search through constrained reaction branches

Working draft, version 0.1 — 13 November 2026. This document accompanies [the main text](MLPATH_MainText.md). It describes saved development cohorts, their provenance, and the scope of their validation. It does not adopt the proposed changes in the 12 November general-rules options document.

## S1. Evidence sources and table reconstruction

Paths in the following table are relative to the repository root. Historical files are inputs and are left unchanged. Earlier versions of the HF-6 summaries are excluded because the b05 symmetry matching was subsequently corrected.

| Evidence | Source artifact | Use |
| --- | --- | --- |
| Local HF/Dimer comparison | `benchmarks/explore_d2/results/2026-11-10/hf6_local_summary_v3.json` | Per-item search status, cost, saddle order, and reference matching |
| End-to-end sA | `benchmarks/explore_d2/results/2026-11-11/hf6_e2e_summary_sA_v3.json` | Recorded budgets and first strict-reference-hit costs |
| End-to-end sB | `benchmarks/explore_d2/results/2026-11-11/hf6_e2e_summary_sB_v2.json` | Second-seed comparison, including the unresolved b05 order |
| Mirror audit | `benchmarks/explore_d2/results/2026-11-12/hf6_mirror_match_check.json` | One additional local N1 Dimer reference match |
| Dimer identity audit | `benchmarks/explore_d2/results/2026-11-12/dimer_identity/summary.json` | Representative grouping, order recheck, and endpoint bond-graph categories |
| Cross-system interpretation | `docs/HF5_SUMMARY_2026-11-10.md` | Sensitivities, noise definitions, and the original N1 connection assignments |
| Corrections and execution records | `docs/DONE_LOG.md`, entries dated 8–12 November | Historical settings, subsequent corrections, execution provenance, and endpoint-assignment discrepancy |

Run the following from the repository root to regenerate the numerical evidence with standard-library Python:

```sh
python paper/tools/build_draft_evidence.py
```

The `--figures` option additionally requires matplotlib. The existing scientific environment supplies it:

```sh
MPLCONFIGDIR=/tmp/mlpath-draft-mpl /home/ruigengji/miniforge3/envs/openmm_dev/bin/python paper/tools/build_draft_evidence.py --figures
```

The script reads five explicit summary files. It does not load a checkpoint, evaluate forces, reconstruct the large journals, reclassify a saddle, or alter a search. The output files are:

| Derived artifact under `paper/` | Contents |
| --- | --- |
| `data/local_results.csv` | Six system/arm aggregate records used by Table 1 |
| `data/local_items.csv` | All 618 item/arm records, with original and draft reference matches |
| `data/end_to_end_results.csv` | Twelve system/seed/arm records used by Table 2 |
| `data/reference_hits.csv` | Thirty-two per-reference observations, including budget-limited non-hits |
| `data/draft_evidence.json` | Source paths, scoring overlay, aggregate evidence, and identity categories |
| `figures/Fig1_local.{png,pdf}` | Total local search cost and strict reference hits |
| `figures/Fig2_recall.{png,pdf}` | Per-seed strict recall versus normalized recorded budget |

## S2. Model, inputs, and numerical settings

### S2.1. Calculator and electronic state

The calculator factory is `benchmarks/b04_sn2_chloride/calculator.py:polar1`. It selects the local file `/home/ruigengji/MLP/mace/MACE-POLAR-1-M.model` with `model_type="PolarMACE"`. The benchmark records contain a historical checkpoint digest, but this draft has not rehashed the model or established its exact correspondence to a public release. The model-family publication cited in the main text does not by itself resolve that correspondence.

All HF-6 searches and the saved validation summaries use float32. b04 and b05 have total charge −1 and multiplicity 1. N1 has total charge 0 and multiplicity 1. Molecular inputs and reference artifacts retain atom identities so that permutations, symmetry, and forming/breaking bonds can be assessed separately.

### S2.2. Local-cohort origins

| System | Frozen HS source under `benchmarks/explore_d2/results/` | D2 rule version | HS inputs |
| --- | --- | --- | ---: |
| b04 | `2026-11-08/b04_hf5_v12_sA/` | D2-HF-v1.2 | 23 |
| b05 | `2026-11-08/b05_hf_v14_sA/` | D2-HF-v1.4 | 140 |
| N1 | `2026-11-09/n1_hf5_v12_sA_n0020/` | D2-HF-v1.2 | 146 |

These are different development configurations. The b05 run enables the later bracket configuration and parallel dispatch; the b04 and N1 local input cohorts retain the earlier setup. All local HS records, including failed items, enter the comparison. Reference structures do not supply search directions.

The end-to-end b04 comparison uses the later b04 v1.4 runs, whereas the b04 local input cohort is v1.2. The sA/sB end-to-end budget and recall values are taken directly from their specified summaries. Thus the local and end-to-end panels answer related but distinct cohort-specific questions.

### S2.3. Shared candidate tolerances and system-specific parameters

| Quantity | Value or policy |
| --- | --- |
| Maximum atomic full-force norm | 0.01 eV Å⁻¹ |
| Candidate transverse-force norm | 0.003 eV Å⁻¹ |
| Fine midpoint transverse-force norm | 0.001 eV Å⁻¹ where configured |
| Coordinate residual | ≤ 10⁻⁸ in the coordinate's units |
| Standard bracket sign-resolution tolerance | 0.003 eV Å⁻¹ for distance coordinates |
| Paired continuation derivative tolerance | 3√2 σ when a measured σ is supplied |
| Local per-arm/per-item budget | 3,000 forward evaluations |
| Reference-match RMSD | ≤ 0.15 Å |
| Reference-match energy difference | ≤ 0.05 eV |
| Post-search analytic-Hessian cross-check tolerance | 0.1 in the shared validation-script configuration |

| Parameter | b04 | b05 | N1 |
| --- | ---: | ---: | ---: |
| Reaction-force σ, eV Å⁻¹ | 0.0042 | 0.0025 | 0.0027 |
| Geometry noise, Å | 0.00062 | 0.007 | 0.0104 |
| Scan energy resolution, eV, principal cohorts | 0.001953125 | 0.00390625 | 0.001953125 |
| Coarse midpoint sign bound, eV Å⁻¹ | 0.0237 in later b04 bracket runs | 0.025 | Not set |
| Coarse midpoint transverse tolerance, eV Å⁻¹ | 0.002 when the coarse configuration is enabled | 0.002 | Not enabled |

The b04 v1.2 local cohort has no coarse midpoint bracket configuration. The later b04 v1.4 runs use the measured b04 sign bound. These values describe executed cohorts; they are not a newly calibrated universal parameter set. The historical noise-measurement procedures differ across systems. In particular, a larger N1 estimate of 0.0527 eV included changes of conformational basin, and an additional run using that estimate is discussed in HF-5 but is not pooled into the main HF-6 local comparison.

### S2.4. Dimer initialization

Local Dimer uses the peak scan geometry within the input bracket and an initial direction from the normalized intent-coordinate gradient. End-to-end DB uses the saved L0 poses, ASE's Gaussian displacement with standard deviation 0.1 Å, and random initial modes under the benchmark's deterministic seed schedule. Both use ASE defaults for the remaining Dimer controls, including the recorded 10⁻⁴ Å separation. The benchmark scripts and saved execution environment govern the historical settings; current upstream defaults are not substituted for them.

The publication and official documentation for the implementation were checked during drafting. The historical runtime package versions still require assembly into the release record. No Dimer-separation sensitivity series is available in the present evidence package.

## S3. Search semantics and validation

### S3.1. Meaning of the recorded statuses

| Record | Evidence supplied |
| --- | --- |
| `converged_candidate` | A returned HS geometry meets the force and coordinate residual contract |
| `unresolved:hysteresis` or branch disagreement | The bounded continuation did not establish the required branch agreement |
| `transverse_not_converged` | The constrained transverse relaxation did not reach its requested threshold |
| `walk_relax_unconverged` | A continuation trial did not complete the relaxation contract |
| `connected` in HC | E/F connection tracking returned its configured endpoint evidence |
| Independent `order == 1` | The selected post-search curvature routine classified the candidate as index one |
| `internal_mode_recheck_disagrees` | A numerical consistency check left the order unresolved |
| QM optimization/frequency/IRC records | A separate level of theory supports the indicated geometry, character, or connection |

These statuses answer different questions. In particular, a branch failure does not diagnose absence of a reaction, and HC connectivity does not automatically establish the candidate's Hessian order.

### S3.2. Hessian cross-check and unresolved order

The HF-6 post-search validation scripts changed the analytic-Hessian cross-check tolerance from 0.01 to 0.1. The recorded rationale was to detect unit/layout-scale discrepancies while avoiding rejection at the observed fp32 cross-check scatter, which reached approximately 0.021 in b04. This revision applies to the validation-script configuration. It did not change the PRRS default or search trajectories.

The separate internal-mode consistency check was retained. Some candidates switch between index-one and unresolved classifications across repeated evaluations. All manuscript counts use the specified saved snapshot, with the subsequent mirror audit overlaid only on matching. The 19 b05 HF local geometric matches include three candidates without a resolved index-one classification, leaving 16 strict hits. In end-to-end sB, one of the three b05 geometric matches is similarly excluded from strict recall.

### S3.3. Symmetry and mirror correction

The older matching implementation used the saddle graph's automorphisms alone. It missed equivalent b05 candidates involving different numbered methyl hydrogens. The selected HF-6 snapshots instead use permutations from the saddle graph and the first reactant-pose graph. The latter supplies equivalent-reactant atom permutations absent from a partially formed saddle graph.

A later audit considered reflection in addition to proper rotation. It found one additional local N1 Dimer match, item 8, to exo s-cis: approximately 0.012 Å after reflection and 0.002 eV energy difference. The draft's N1 local Dimer count is therefore three strict attempts spanning two reference classes. End-to-end recall is unchanged. Original summary files and their original matches are preserved, and `local_items.csv` retains both match columns.

Distinct-saddle counts originally obtained with proper-rotation-only grouping have not been recomputed with a comprehensive mirror policy. The main table consequently reports strict reference-hit attempts and avoids treating those grouping counts as a finalized count of distinct channels.

## S4. Cost accounting and observation limits

The local cost is the sum of all recorded forward evaluations across each arm's frozen input set. It includes budget-exhausted and otherwise failed attempts. It excludes construction of the frozen scans and post-search saddle validation. The ratios 42,048/4,360 = 9.64 and 371,593/63,659 = 5.84 are total local-cohort expenditure ratios. They are not wall-time speedups or per-unique-channel efficiencies.

End-to-end DB receives a budget equal to the observed D2 ledger total. This ledger boundary includes D2's charged scan and connection work; DB consumes its budget on Dimer attempts from the shared saved L0 library. The archived data do not yet establish a common complete-job cost from raw reactant files through fully certified endpoints. The post-search validation ledger and later Dimer identity audit are outside the plotted search budgets.

| System | Seed | D2 first strict hit | Recorded budget | First-hit budget fraction |
| --- | --- | ---: | ---: | ---: |
| b04 | sA | 24,800 | 26,129 | 0.9491 |
| b04 | sB | 26,506 | 28,071 | 0.9442 |
| b05 | sA | 348,119 | 408,888 | 0.8514 |
| b05 | sB | 348,388 | 405,075 | 0.8601 |
| N1 | sA | No designated-reference hit | 376,296 | Unobserved |
| N1 | sB | No designated-reference hit | 377,745 | Unobserved |

The b05 sB geometric match at cost 348,076 is unresolved in order and does not establish the strict first-hit time. Individual reference costs are retained in `reference_hits.csv`. Two seeds provide descriptive replication, without a calibrated population success rate or uncertainty interval. The 309 local intents are nested within three systems and must not be treated as 309 independent reaction systems.

## S5. Reference-external candidate identity audit

The saved audit selects Dimer candidates that were index one and nonmatching in the selected snapshots, groups them using the existing geometric/energy rule, and rechecks one representative per group. For index-one representatives, it displaces along both signs of the unstable mode with a maximum atomic displacement of 0.15 Å and uses ASE BFGS relaxation to 0.01 eV Å⁻¹, with a 1,000-step cap. Endpoints are classified by bond graphs. No endpoint Hessian minimum certification or relay extension was performed in this audit.

| System | Selected candidates | Groups | Same endpoint bond graph | One reactant endpoint and different endpoint graphs | Neither endpoint reactant, with different graphs | Order unresolved on recheck | Incomplete descent |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| b04 | 1 | 1 | 1 | 0 | 0 | 0 | 0 |
| b05 | 42 | 24 | 9 | 1 | 11 | 3 | 0 |
| N1 | 37 | 37 | 8 | 10 | 5 | 13 | 1 |

For b05 group 0, the six members comprise four local observations and two DB sB observations. One descending side has the reactant graph; the other breaks C0–Cl1 and forms C0–O14. Its nearest D2-event RMSD is approximately 0.445 Å. The group supports a provisional substitution connection on the model surface, with endpoint minimum and QM validation still absent.

For N1, four audited saddle groups match D2 event geometries. One shared saddle was assigned a two-bond [2+2]-type endpoint in the later descent audit, although the HF-5 HC endpoint record assigned a single-bond change. The procedures yield different endpoint assignments and have not been reconciled. The draft therefore avoids a categorical claim that all D2 N1 events are mechanistically single-bond stepwise reactions. None matches the four concerted DA reference geometries.

The audit cost was recorded as 11,352 forwards and 1,488 s. It pools local and end-to-end observations for identity analysis; those grouped discoveries do not form an additional independently sampled performance cohort.

## S6. QM evidence and model-surface claims

The b04 historical record is `docs/experiments/b04_qm/README.md`. It describes QM refinement from a PRRS-generated POLAR-1 saddle, r2SCAN-3c frequencies and two-sided IRC, endpoint optimization/frequencies, and ωB97M-V/def2-TZVPD single points. This validates the historical reference event. Its provenance remains the earlier PRRS calculation rather than a new independent QM calculation for every HF candidate.

The available b05 summary, `docs/experiments/b05_qm/summary.json`, records optimized TS and endpoint structures for a recovered geometry and a control. Both TS entries report one imaginary frequency and the endpoint entries report none. This draft uses those records as representative QM evidence, without reprocessing the underlying frequency output blocks or asserting a full QM IRC connection for every b05 reference.

N1's four designated DA references have model-surface saddle and two-sided-descent evidence described in the reference-building records. They lack a QM validation axis in this study. The additional b05 substitution group also lacks QM evidence. Relative energies between model saddles are reported as model-surface observations; they do not establish free-energy barriers, kinetics, or selectivity.

## S7. Sensitivity and transferability boundaries

The b04 development runs include repeated execution, changed seed, rigid rotation, and a smaller scan step. All six examined runs recovered the designated SN2 geometry, although some historical order classifications depended on the validation threshold. The count of promoted ridges changed with rotation and scan step.

The same-graph ridge heuristic compares adjacent peak-region bond graphs and reduces local-search allocation to apparently unchanged graphs. b05 supplied the original evidence for this expenditure reduction, but b04 and N1 supplied counterexamples. In a later offline replay, comparing the scan start and the relaxed post-peak low point recovered a more useful partition on five existing runs. This is an analysis of previously observed trajectories. It has not been validated as a new frozen search rule on independent runs.

The three proposed general-rule revisions—ridge classification, fixed numerical calibration, and energy-noise definition—remain options in `docs/HF_GENERAL_RULES_OPTIONS_2026-11-12.md`. They were not applied retroactively to these results. A later manuscript version should retain this development cohort separately from any validation under revised rules.

## S8. Literature provenance and release status

The foundational Dimer and CI-NEB metadata were checked against the University of Iceland author records; FIRE against the APS article page; ASE against its author-repository record and official Dimer documentation; and MACE-POLAR-1 against arXiv:2602.19411v1. These sources support the specific method and model-family descriptions in the main text. This is a starting bibliography, not a completed survey of related branch-following and reaction-search methods.

The local model's release correspondence, license and redistribution conditions, full historical runtime versions, final authorship/funding statements, public dataset deposit, and submission venue remain unassembled. These release items do not alter the provenance of the saved development results, but they are needed for a submission-ready package.
