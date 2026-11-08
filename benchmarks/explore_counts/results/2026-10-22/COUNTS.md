# L0 pure-geometry counts

Rule version: `L0-v1`. Pure geometry, no model loaded. Counting only; no geometry here is a proposal.

Parameters:

- steps_per_turn: 6
- conformer cap: 1000000
- seed: 0
- n_samples: [8, 32]
- shells: measured f of the input, f-0.4, f+0.4 (the +-0.4 are sensitivity only, not proposals)
- intent poses per (system, shell, n_samples=32): 16
- n_far: 0
- pair locality rules: ['shared_atom', 'adjacent', 'all']
- scan points per intent: 18
- forwards per point: 8.7
- seconds per point: 1.07
- systems: ['b04', 'b05', 'n1', 'b09']
- grid0_note: Per-pose intent counts come from grid-0 conformer-pair poses only (the first 16 placed poses); est. poses/intents = placed poses x conformer pairs, i.e. they assume every conformer pair has the same per-pose count (matters for b09 x210, n1 x6)
- pair_cost_assumption: ASSUMPTION: every pair intent is charged the same 18 points x 1.07 s as a single-distance scan; a pair intent may cost more

Cost note: b04 (6 atoms) measured, extrapolated; larger systems likely slower per point (more BFGS steps, slower forward); order of magnitude only.

Note: Per-pose intent counts come from grid-0 conformer-pair poses only (the first 16 placed poses); est. poses/intents = placed poses x conformer pairs, i.e. they assume every conformer pair has the same per-pose count (matters for b09 x210, n1 x6).

Note: ASSUMPTION: every pair intent is charged the same 18 points x 1.07 s as a single-distance scan; a pair intent may cost more.

## b04

- atoms 6, charge -1, fragments 2
- measured contact ratio f of the input: 2.0225

Conformers (steps_per_turn 6, cap 1e6)

| fragment | genuine torsions | rotors not enumerated | grid | generated | status counts |
|---|---|---|---|---|---|
| CH3Cl|2b70d83c1c80f735,3e9540dc663201cf,c4e3929869b32963,c4e3929869b32963,c4e3929869b32963 | 0 | 0 | 1 | 1 | {'generated': 1} |
| Cl|e1734f3e5568970e | 0 | 0 | 1 | 1 | {'generated': 1} |

Conformer pairs (generated x generated): 1

| shell | n_samples | poses | placed | status | singles/pose (bonded, near; intra, inter) [grid-0 poses] | pair intents/pose (shared_atom, adjacent, all) [grid-0 poses] | far classes/pose [grid-0 poses] | est. poses (placed x conformer pairs) | est. intents single (assumes equal per conformer pair) | single+pair adjacent (same assumption) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.622 (f-0.4 (sensitivity only, not a proposal)) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 1.622 (f-0.4 (sensitivity only, not a proposal)) | 32 | 32 | 32 | {'placed': 32} | 14.6 (8.0, 6.6; 8.0, 6.6) | 66, 98, 98 | 1.1 | 32 | 468 | 3,588 |
| 2.022 (measured f) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 2.022 (measured f) | 32 | 32 | 32 | {'placed': 32} | 11.1 (8.0, 3.1; 8.0, 3.1) | 42, 55, 55 | 2.4 | 32 | 356 | 2,124 |
| 2.422 (f+0.4 (sensitivity only, not a proposal)) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 2.422 (f+0.4 (sensitivity only, not a proposal)) | 32 | 32 | 32 | {'placed': 32} | 9.0 (8.0, 1.0; 8.0, 1.0) | 26, 32, 32 | 3.0 | 32 | 288 | 1,320 |

Cost extrapolation (b04 (6 atoms) measured, extrapolated; larger systems likely slower per point (more BFGS steps, slower forward); order of magnitude only)

| shell | n_samples | rule | intents | forwards | wall hours | GPU-hours (1x 2080 Ti) |
|---|---|---|---|---|---|---|
| f-0.4 (sensitivity only, not a proposal) 1.622 | 32 | single | 468 | 7.329e+04 | 2.50 | 2.50 |
| f-0.4 (sensitivity only, not a proposal) 1.622 | 32 | single+pair[shared_atom] | 2,580 | 4.040e+05 | 13.80 | 13.80 |
| f-0.4 (sensitivity only, not a proposal) 1.622 | 32 | single+pair[adjacent] | 3,588 | 5.619e+05 | 19.20 | 19.20 |
| f-0.4 (sensitivity only, not a proposal) 1.622 | 32 | single+pair[all] | 3,588 | 5.619e+05 | 19.20 | 19.20 |
| measured f 2.022 | 32 | single | 356 | 5.575e+04 | 1.90 | 1.90 |
| measured f 2.022 | 32 | single+pair[shared_atom] | 1,716 | 2.687e+05 | 9.18 | 9.18 |
| measured f 2.022 | 32 | single+pair[adjacent] | 2,124 | 3.326e+05 | 11.36 | 11.36 |
| measured f 2.022 | 32 | single+pair[all] | 2,124 | 3.326e+05 | 11.36 | 11.36 |
| f+0.4 (sensitivity only, not a proposal) 2.422 | 32 | single | 288 | 4.510e+04 | 1.54 | 1.54 |
| f+0.4 (sensitivity only, not a proposal) 2.422 | 32 | single+pair[shared_atom] | 1,128 | 1.766e+05 | 6.03 | 6.03 |
| f+0.4 (sensitivity only, not a proposal) 2.422 | 32 | single+pair[adjacent] | 1,320 | 2.067e+05 | 7.06 | 7.06 |
| f+0.4 (sensitivity only, not a proposal) 2.422 | 32 | single+pair[all] | 1,320 | 2.067e+05 | 7.06 | 7.06 |
ASSUMPTION: every pair intent is charged the same 18 points x 1.07 s as a single-distance scan; a pair intent may cost more.

Per-pose intent counts come from grid-0 conformer-pair poses only (the first 16 placed poses); est. poses/intents = placed poses x conformer pairs, i.e. they assume every conformer pair has the same per-pose count (matters for b09 x210, n1 x6).

## b05

- atoms 16, charge -1, fragments 2
- measured contact ratio f of the input: 2.4008
- LOWER BOUND: conformer counts exclude non-enumerated rotors 3 in C4H9Cl; true conformer counts are larger (ring conformations are never enumerated in v1)

Conformers (steps_per_turn 6, cap 1e6)

| fragment | genuine torsions | rotors not enumerated | grid | generated | status counts |
|---|---|---|---|---|---|
| C4H9Cl|69bec734e5749d86,7b6fcaad228d7a4c,86d490de8c338b11,86d490de8c338b11,86d490de8c338b11,86d490de8c338b11,86d490de8c338b11,86d490de8c338b11,86d490de8c338b11,86d490de8c338b11,86d490de8c338b11,fec7f36755270831,fec7f36755270831,fec7f36755270831 | 0 | 3 | 1 | 1 | {'generated': 1} |
| HO|20e0501a2ce64366,e7ef864bc9081975 | 0 | 0 | 1 | 1 | {'generated': 1} |

Conformer pairs (generated x generated): 1

| shell | n_samples | poses | placed | status | singles/pose (bonded, near; intra, inter) [grid-0 poses] | pair intents/pose (shared_atom, adjacent, all) [grid-0 poses] | far classes/pose [grid-0 poses] | est. poses (placed x conformer pairs) | est. intents single (assumes equal per conformer pair) | single+pair adjacent (same assumption) |
|---|---|---|---|---|---|---|---|---|---|---|
| 2.001 (f-0.4 (sensitivity only, not a proposal)) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 2.001 (f-0.4 (sensitivity only, not a proposal)) | 32 | 32 | 32 | {'placed': 32} | 97.9 (28.0, 69.9; 76.0, 21.9) | 1,136, 2,270, 4,768 | 6.8 | 32 | 3,132 | 75,772 |
| 2.401 (measured f) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 2.401 (measured f) | 32 | 32 | 32 | {'placed': 32} | 90.5 (28.0, 62.5; 76.0, 14.5) | 972, 1,933, 4,050 | 8.4 | 32 | 2,896 | 64,760 |
| 2.801 (f+0.4 (sensitivity only, not a proposal)) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 2.801 (f+0.4 (sensitivity only, not a proposal)) | 32 | 32 | 32 | {'placed': 32} | 84.5 (28.0, 56.5; 76.0, 8.5) | 850, 1,674, 3,505 | 8.9 | 32 | 2,704 | 56,256 |

Cost extrapolation (b04 (6 atoms) measured, extrapolated; larger systems likely slower per point (more BFGS steps, slower forward); order of magnitude only)

| shell | n_samples | rule | intents | forwards | wall hours | GPU-hours (1x 2080 Ti) |
|---|---|---|---|---|---|---|
| f-0.4 (sensitivity only, not a proposal) 2.001 | 32 | single | 3,132 | 4.905e+05 | 16.76 | 16.76 |
| f-0.4 (sensitivity only, not a proposal) 2.001 | 32 | single+pair[shared_atom] | 39,468 | 6.181e+06 | 211.15 | 211.15 |
| f-0.4 (sensitivity only, not a proposal) 2.001 | 32 | single+pair[adjacent] | 75,772 | 1.187e+07 | 405.38 | 405.38 |
| f-0.4 (sensitivity only, not a proposal) 2.001 | 32 | single+pair[all] | 155,724 | 2.439e+07 | 833.12 | 833.12 |
| measured f 2.401 | 32 | single | 2,896 | 4.535e+05 | 15.49 | 15.49 |
| measured f 2.401 | 32 | single+pair[shared_atom] | 33,992 | 5.323e+06 | 181.86 | 181.86 |
| measured f 2.401 | 32 | single+pair[adjacent] | 64,760 | 1.014e+07 | 346.47 | 346.47 |
| measured f 2.401 | 32 | single+pair[all] | 132,512 | 2.075e+07 | 708.94 | 708.94 |
| f+0.4 (sensitivity only, not a proposal) 2.801 | 32 | single | 2,704 | 4.234e+05 | 14.47 | 14.47 |
| f+0.4 (sensitivity only, not a proposal) 2.801 | 32 | single+pair[shared_atom] | 29,920 | 4.685e+06 | 160.07 | 160.07 |
| f+0.4 (sensitivity only, not a proposal) 2.801 | 32 | single+pair[adjacent] | 56,256 | 8.810e+06 | 300.97 | 300.97 |
| f+0.4 (sensitivity only, not a proposal) 2.801 | 32 | single+pair[all] | 114,864 | 1.799e+07 | 614.52 | 614.52 |
ASSUMPTION: every pair intent is charged the same 18 points x 1.07 s as a single-distance scan; a pair intent may cost more.

Per-pose intent counts come from grid-0 conformer-pair poses only (the first 16 placed poses); est. poses/intents = placed poses x conformer pairs, i.e. they assume every conformer pair has the same per-pose count (matters for b09 x210, n1 x6).

## n1

- atoms 19, charge 0, fragments 2
- NOTE: MMFF geometry, counting only, not a POLAR-1 minimum; the contact ratio is an artefact of placing the MMFF fragments 8 A apart, not a complex contact ratio
- placement f (8 A artificial, not a complex contact ratio): 3.6674
- LOWER BOUND: conformer counts exclude non-enumerated rotors 1 in C3H4O; true conformer counts are larger (ring conformations are never enumerated in v1)

Conformers (steps_per_turn 6, cap 1e6)

| fragment | genuine torsions | rotors not enumerated | grid | generated | status counts |
|---|---|---|---|---|---|
| C3H4O|045321b2a7074753,527e2029728d6f31,8e188e7ac74b0248,8fbe24ba235ac7a7,9a4d550cdcff4f84,9a4d550cdcff4f84,b2e0df8e05daf177,f661a844ddab4f7a | 1 | 1 | 6 | 6 | {'generated': 6} |
| C5H6|00e24f2b7b923640,00e24f2b7b923640,57e895b0f4df58ba,57e895b0f4df58ba,6058023ca7356ad0,6058023ca7356ad0,632ff7d0ab4f4d75,632ff7d0ab4f4d75,e3f2d767ccd57a61,f9994bdb4f4b194d,f9994bdb4f4b194d | 0 | 0 | 1 | 1 | {'generated': 1} |

Conformer pairs (generated x generated): 6

| shell | n_samples | poses | placed | status | singles/pose (bonded, near; intra, inter) [grid-0 poses] | pair intents/pose (shared_atom, adjacent, all) [grid-0 poses] | far classes/pose [grid-0 poses] | est. poses (placed x conformer pairs) | est. intents single (assumes equal per conformer pair) | single+pair adjacent (same assumption) |
|---|---|---|---|---|---|---|---|---|---|---|
| 3.267 (f-0.4 (sensitivity only, not a proposal)) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 3.267 (f-0.4 (sensitivity only, not a proposal)) | 32 | 32 | 32 | {'placed': 32} | 45.5 (36.0, 9.5; 42.0, 3.5) | 169, 378, 992 | 45.8 | 192 | 8,736 | 81,408 |
| 3.667 (placement f (8 A artificial)) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 3.667 (placement f (8 A artificial)) | 32 | 32 | 32 | {'placed': 32} | 42.6 (36.0, 6.6; 42.0, 0.6) | 145, 327, 867 | 46.0 | 192 | 8,184 | 70,920 |
| 4.067 (f+0.4 (sensitivity only, not a proposal)) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 4.067 (f+0.4 (sensitivity only, not a proposal)) | 32 | 32 | 32 | {'placed': 32} | 42.0 (36.0, 6.0; 42.0, 0.0) | 140, 316, 840 | 46.0 | 192 | 8,064 | 68,736 |

Cost extrapolation (b04 (6 atoms) measured, extrapolated; larger systems likely slower per point (more BFGS steps, slower forward); order of magnitude only)

| shell | n_samples | rule | intents | forwards | wall hours | GPU-hours (1x 2080 Ti) |
|---|---|---|---|---|---|---|
| f-0.4 (sensitivity only, not a proposal) 3.267 | 32 | single | 8,736 | 1.368e+06 | 46.74 | 46.74 |
| f-0.4 (sensitivity only, not a proposal) 3.267 | 32 | single+pair[shared_atom] | 41,232 | 6.457e+06 | 220.59 | 220.59 |
| f-0.4 (sensitivity only, not a proposal) 3.267 | 32 | single+pair[adjacent] | 81,408 | 1.275e+07 | 435.53 | 435.53 |
| f-0.4 (sensitivity only, not a proposal) 3.267 | 32 | single+pair[all] | 199,296 | 3.121e+07 | 1066.23 | 1066.23 |
| placement f (8 A artificial) 3.667 | 32 | single | 8,184 | 1.282e+06 | 43.78 | 43.78 |
| placement f (8 A artificial) 3.667 | 32 | single+pair[shared_atom] | 36,024 | 5.641e+06 | 192.73 | 192.73 |
| placement f (8 A artificial) 3.667 | 32 | single+pair[adjacent] | 70,920 | 1.111e+07 | 379.42 | 379.42 |
| placement f (8 A artificial) 3.667 | 32 | single+pair[all] | 174,600 | 2.734e+07 | 934.11 | 934.11 |
| f+0.4 (sensitivity only, not a proposal) 4.067 | 32 | single | 8,064 | 1.263e+06 | 43.14 | 43.14 |
| f+0.4 (sensitivity only, not a proposal) 4.067 | 32 | single+pair[shared_atom] | 34,944 | 5.472e+06 | 186.95 | 186.95 |
| f+0.4 (sensitivity only, not a proposal) 4.067 | 32 | single+pair[adjacent] | 68,736 | 1.076e+07 | 367.74 | 367.74 |
| f+0.4 (sensitivity only, not a proposal) 4.067 | 32 | single+pair[all] | 169,344 | 2.652e+07 | 905.99 | 905.99 |
ASSUMPTION: every pair intent is charged the same 18 points x 1.07 s as a single-distance scan; a pair intent may cost more.

Per-pose intent counts come from grid-0 conformer-pair poses only (the first 16 placed poses); est. poses/intents = placed poses x conformer pairs, i.e. they assume every conformer pair has the same per-pose count (matters for b09 x210, n1 x6).

## b09

- atoms 29, charge 0, fragments 2
- NOTE: biased toward TS-2OMe; counting only, not a blind-search start
- measured contact ratio f of the input: 2.0551
- LOWER BOUND: conformer counts exclude non-enumerated rotors 2 in C3H6O2; true conformer counts are larger (ring conformations are never enumerated in v1)

Conformers (steps_per_turn 6, cap 1e6)

| fragment | genuine torsions | rotors not enumerated | grid | generated | status counts |
|---|---|---|---|---|---|
| C3H6O2|160add52eb764b0e,160add52eb764b0e,160add52eb764b0e,160add52eb764b0e,160add52eb764b0e,160add52eb764b0e,8a417c8fb2244bc8,a405b93c64b249c6,a405b93c64b249c6,cb8466acff5d1b7d,cb8466acff5d1b7d | 2 | 2 | 36 | 35 | {'generated': 35, 'graph_changed': 1} |
| C7H10O|01256f6732f4c71d,0354ffd8f06c9f07,1052608a33bf4096,1f924b524573e6f8,1f924b524573e6f8,294fb5ecf1c18af5,2bd822de086ae05b,2bd822de086ae05b,55a8b9e06d0ade61,5d2fbfb851416e6c,5d2fbfb851416e6c,761cf636d30cacc4,813ab01d5d1f0cc6,845b96621157d7a9,983c5bcc9ceb29c6,999d760494cdcc5e,e659b45834915507,e659b45834915507 | 1 | 0 | 6 | 6 | {'generated': 6} |

Conformer pairs (generated x generated): 210

| shell | n_samples | poses | placed | status | singles/pose (bonded, near; intra, inter) [grid-0 poses] | pair intents/pose (shared_atom, adjacent, all) [grid-0 poses] | far classes/pose [grid-0 poses] | est. poses (placed x conformer pairs) | est. intents single (assumes equal per conformer pair) | single+pair adjacent (same assumption) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.655 (f-0.4 (sensitivity only, not a proposal)) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 1.655 (f-0.4 (sensitivity only, not a proposal)) | 32 | 32 | 32 | {'placed': 32} | 216.2 (56.0, 160.2; 126.0, 90.2) | 3,488, 8,893, 23,764 | 75.6 | 6,720 | 1,453,200 | 61,214,160 |
| 2.055 (measured f) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 2.055 (measured f) | 32 | 32 | 32 | {'placed': 32} | 185.2 (56.0, 129.2; 126.0, 59.2) | 2,486, 6,136, 17,298 | 78.4 | 6,720 | 1,244,880 | 42,477,120 |
| 2.455 (f+0.4 (sensitivity only, not a proposal)) | 8 | 8 | 8 | {'placed': 8} | n/a | n/a | n/a | n/a | n/a | n/a |
| 2.455 (f+0.4 (sensitivity only, not a proposal)) | 32 | 32 | 32 | {'placed': 32} | 159.6 (56.0, 103.6; 126.0, 33.6) | 1,758, 4,242, 12,710 | 80.0 | 6,720 | 1,072,680 | 29,578,920 |

Cost extrapolation (b04 (6 atoms) measured, extrapolated; larger systems likely slower per point (more BFGS steps, slower forward); order of magnitude only)

| shell | n_samples | rule | intents | forwards | wall hours | GPU-hours (1x 2080 Ti) |
|---|---|---|---|---|---|---|
| f-0.4 (sensitivity only, not a proposal) 1.655 | 32 | single | 1,453,200 | 2.276e+08 | 7774.62 | 7774.62 |
| f-0.4 (sensitivity only, not a proposal) 1.655 | 32 | single+pair[shared_atom] | 24,892,560 | 3.898e+09 | 133175.20 | 133175.20 |
| f-0.4 (sensitivity only, not a proposal) 1.655 | 32 | single+pair[adjacent] | 61,214,160 | 9.586e+09 | 327495.76 | 327495.76 |
| f-0.4 (sensitivity only, not a proposal) 1.655 | 32 | single+pair[all] | 161,150,640 | 2.524e+10 | 862155.92 | 862155.92 |
| measured f 2.055 | 32 | single | 1,244,880 | 1.949e+08 | 6660.11 | 6660.11 |
| measured f 2.055 | 32 | single+pair[shared_atom] | 17,947,440 | 2.811e+09 | 96018.80 | 96018.80 |
| measured f 2.055 | 32 | single+pair[adjacent] | 42,477,120 | 6.652e+09 | 227252.59 | 227252.59 |
| measured f 2.055 | 32 | single+pair[all] | 117,487,440 | 1.840e+10 | 628557.80 | 628557.80 |
| f+0.4 (sensitivity only, not a proposal) 2.455 | 32 | single | 1,072,680 | 1.680e+08 | 5738.84 | 5738.84 |
| f+0.4 (sensitivity only, not a proposal) 2.455 | 32 | single+pair[shared_atom] | 12,883,080 | 2.017e+09 | 68924.48 | 68924.48 |
| f+0.4 (sensitivity only, not a proposal) 2.455 | 32 | single+pair[adjacent] | 29,578,920 | 4.632e+09 | 158247.22 | 158247.22 |
| f+0.4 (sensitivity only, not a proposal) 2.455 | 32 | single+pair[all] | 86,480,520 | 1.354e+10 | 462670.78 | 462670.78 |
ASSUMPTION: every pair intent is charged the same 18 points x 1.07 s as a single-distance scan; a pair intent may cost more.

Per-pose intent counts come from grid-0 conformer-pair poses only (the first 16 placed poses); est. poses/intents = placed poses x conformer pairs, i.e. they assume every conformer pair has the same per-pose count (matters for b09 x210, n1 x6).

## Summary (n_samples 32, measured-f and sensitivity shells; per-pose intents from first 16 placed poses)

| system | note | shell | placed poses | singles/pose | single+pair(adjacent)/pose | est. intents single | est. intents single+adjacent | GPU-h single | GPU-h single+adjacent |
|---|---|---|---|---|---|---|---|---|---|
| b04 | - | 1.622 (f-0.4, sensitivity only) | 32 | 14.6 | 112.1 | 468 | 3,588 | 2.50 | 19.20 |
| b04 | - | 2.022 (measured f) | 32 | 11.1 | 66.4 | 356 | 2,124 | 1.90 | 11.36 |
| b04 | - | 2.422 (f+0.4, sensitivity only) | 32 | 9.0 | 41.2 | 288 | 1,320 | 1.54 | 7.06 |
| b05 | -; conformer counts are lower bounds (rotors not enumerated) | 2.001 (f-0.4, sensitivity only) | 32 | 97.9 | 2,367.9 | 3,132 | 75,772 | 16.76 | 405.38 |
| b05 | -; conformer counts are lower bounds (rotors not enumerated) | 2.401 (measured f) | 32 | 90.5 | 2,023.8 | 2,896 | 64,760 | 15.49 | 346.47 |
| b05 | -; conformer counts are lower bounds (rotors not enumerated) | 2.801 (f+0.4, sensitivity only) | 32 | 84.5 | 1,758.0 | 2,704 | 56,256 | 14.47 | 300.97 |
| n1 | MMFF geometry, counting only, not a POLAR-1 minimum; the contact ratio is an artefact of placing the MMFF fragments 8 A apart, not a complex contact ratio; conformer counts are lower bounds (rotors not enumerated) | 3.267 (f-0.4, sensitivity only) | 32 | 45.5 | 424.0 | 8,736 | 81,408 | 46.74 | 435.53 |
| n1 | MMFF geometry, counting only, not a POLAR-1 minimum; the contact ratio is an artefact of placing the MMFF fragments 8 A apart, not a complex contact ratio; conformer counts are lower bounds (rotors not enumerated) | 3.667 (placement f (8 A artificial)) | 32 | 42.6 | 369.4 | 8,184 | 70,920 | 43.78 | 379.42 |
| n1 | MMFF geometry, counting only, not a POLAR-1 minimum; the contact ratio is an artefact of placing the MMFF fragments 8 A apart, not a complex contact ratio; conformer counts are lower bounds (rotors not enumerated) | 4.067 (f+0.4, sensitivity only) | 32 | 42.0 | 358.0 | 8,064 | 68,736 | 43.14 | 367.74 |
| b09 | biased toward TS-2OMe; counting only, not a blind-search start; conformer counts are lower bounds (rotors not enumerated) | 1.655 (f-0.4, sensitivity only) | 32 | 216.2 | 9,109.2 | 1,453,200 | 61,214,160 | 7774.62 | 327495.76 |
| b09 | biased toward TS-2OMe; counting only, not a blind-search start; conformer counts are lower bounds (rotors not enumerated) | 2.055 (measured f) | 32 | 185.2 | 6,321.0 | 1,244,880 | 42,477,120 | 6660.11 | 227252.59 |
| b09 | biased toward TS-2OMe; counting only, not a blind-search start; conformer counts are lower bounds (rotors not enumerated) | 2.455 (f+0.4, sensitivity only) | 32 | 159.6 | 4,401.6 | 1,072,680 | 29,578,920 | 5738.84 | 158247.22 |

Cost note: b04 (6 atoms) measured, extrapolated; larger systems likely slower per point (more BFGS steps, slower forward); order of magnitude only.

ASSUMPTION: every pair intent is charged the same 18 points x 1.07 s as a single-distance scan; a pair intent may cost more.

Per-pose intent counts come from grid-0 conformer-pair poses only (the first 16 placed poses); est. poses/intents = placed poses x conformer pairs, i.e. they assume every conformer pair has the same per-pose count (matters for b09 x210, n1 x6).
