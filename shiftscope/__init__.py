"""ShiftScope: compare two single-cell populations and localize the shift.

Import the modules you need, e.g. `from shiftscope import io, embed, compare`.
The pipeline order is: io -> embed -> compare -> localize -> drivers -> interpret -> app.

`prioritize` is an analysis branch on top of `compare.rank` (or any hit table): it ranks
screen hits by strong-phenotype x under-studied (grounded in real PubMed counts + a Claude
verdict) into a "validate these" shortlist. `calibration` stress-tests the E-test's detection
limits (power vs cell number and vs effect size) on a known shift.
"""
