"""Mechanical evaluation of the frozen b04 criteria. NOT WRITTEN YET.

Deliberately a stub. Writing checks before the criteria are frozen inverts the order
this project runs on: docs/experiments/p2_check.py opens with "No criterion is computed
here that was not written down in docs/P2_OXOBUTANAL_ACCEPTANCE.md first", and that
sentence is the whole discipline. Code written first would become the criteria by
default, and nobody would notice that they had never been frozen.

When CRITERIA.md is frozen, model this file on docs/experiments/p2_check.py: quote each
criterion's wording next to its check, so a drift between the document and the code is
visible rather than silent.

    python benchmarks/b04_sn2_chloride/check.py runs/<the run>/network.json
"""

import sys

print(__doc__.strip())
print()
print("b04 criteria are not frozen. See CRITERIA.md for what has to be true first,")
print("and preflight.md for the two gates that are currently red.")
sys.exit(2)
