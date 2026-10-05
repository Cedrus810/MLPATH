# Fixed geometry fixtures

These five small extxyz files are byte-for-byte copies of the historical
geometries below, captured on 2026-10-04 for offline reaction-identity and
proposal-priority regressions. Tests must not require the ignored `runs/`
directory for these fixed inputs. Energies/forces in the historical files
are retained as provenance; this is not a new model evaluation.

| Fixture | Original path |
| --- | --- |
| e2e_ethanol/input.extxyz | runs/e2e_ethanol/input.extxyz |
| p1_malonaldehyde/input.extxyz | runs/p1_malonaldehyde/input.extxyz |
| p1_malonaldehyde/c0000_m0000.extxyz | runs/p1_malonaldehyde/c0000_m0000.extxyz |
| p1_malonaldehyde/c0001_m0000.extxyz | runs/p1_malonaldehyde/c0001_m0000.extxyz |
| p2_oxobutanal/p2_input.extxyz | runs/p2_oxobutanal/p2_input.extxyz |

Corpus-level schema/digest gates still require the historical runs. Only
absence of the entire runs directory skips the schema ledger tests;
an existing incomplete corpus continues to fail. Supply the original
corpus to reproduce the recorded golden-sample validation.
