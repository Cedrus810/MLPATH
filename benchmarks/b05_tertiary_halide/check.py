"""Mechanical evaluation of the frozen b05 criteria. NOT WRITTEN YET.

Deliberately a stub, for the same reason b04's is: writing checks before the criteria
are frozen inverts the order this project runs on. docs/experiments/p2_check.py opens
with "No criterion is computed here that was not written down in
docs/P2_OXOBUTANAL_ACCEPTANCE.md first", and that sentence is the whole discipline.
Code written first becomes the criteria by default, and nobody notices they were never
frozen.

When CRITERIA.md is frozen, model this file on docs/experiments/p2_check.py: quote each
criterion's wording next to its check, so drift between the document and the code is
visible rather than silent. This grid has one requirement b04 does not: every "not
found" verdict has to cite that run's own proposal budget, because "never proposed" and
"searched and not found" are different results (../README.md rule 4) and preflight
already found one coordinate at risk of the confusion (r_C_H_beta, gate 6).

    python benchmarks/b05_tertiary_halide/check.py runs/<the run>/network.json
"""

import sys

print(__doc__.strip())
print()
print("b05 criteria are not frozen. See CRITERIA.md for what has to be true first;")
print("prerequisite 2 there is b04's own criteria, which are also not frozen.")
sys.exit(2)
