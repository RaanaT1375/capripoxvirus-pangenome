# Stage 29 — reconciling the absent temporal signal with the published genus rate

This is Option B from `novelty_options.md`, and it turned out to be a
correction rather than an addition.

## The error it fixes

The Results state that at mu = 1.0e-5 substitutions/site/year the 110,321-bp
core alignment is expected to accumulate "about 5 substitutions over 25 years
and about 15 over the 71-year span".

    1.0e-5 x 110,321 x 25 = 27.6
    1.0e-5 x 110,321 x 71 = 78.3

The quoted values are low by a factor of 5.5, and the error runs in the
direction that supports the manuscript's argument. The two figures are
internally consistent with a rate of 1.8e-6, which is 5.5x below the rate the
Methods cite, so one rate was applied consistently — just not the stated one.

Consequence for the argument: structural divergence across Node_147 is
100-250 substitutions, so the excess over the 25-year temporal expectation
falls from the implied 20-50x to 4-9x. Structure still dominates time; the
margin is a few-fold, not an order of magnitude, and the text must say so.

## What it adds

"No temporal signal" and "no detectable temporal signal" are different
claims. Only the second is supported unless the power to detect a slope of mu
is computed, and a reviewer is entitled to ask. The script computes:

- expected substitutions over any span, across a range of rates
- the root-to-tip regression refitted from the stored per-isolate data, which
  also verifies that the table is the one the manuscript used
- a date-randomisation test, the accepted way to report absent temporal signal
- the minimum slope detectable at 80% power, and the power to detect mu
- the standing diversity converted to a timescale (pi / 2mu)

## Reading the power result

Low power (< 50%) is the safe outcome: the absence of signal is consistent
with the published rate, and the finding becomes "these data cannot resolve a
clock of this strength" rather than "no clock exists".

High power is the awkward outcome: a clock at the published rate should have
been visible and is not, which needs an explanation — recombination
disrupting root-to-tip additivity, imprecise collection dates, or a
within-LSDV rate below the genus estimate. It must not then be reported as
unsurprising.

Either way the power figure is an UPPER bound, because root-to-tip distances
share branches and are not independent, so the effective sample size is below
n. A low-power verdict is therefore robust and a high-power verdict is
provisional.

## Run

    conda activate iqtree_env
    cd $PANG/revision_2026-11/29_Temporal_Reconciliation
    python3 scripts/01_temporal_power.py

If it cannot find the per-isolate temporal table it prints a search command;
re-run with TEMPORAL=/path/to/file.csv. Parameters MU, ALN, PI, EPOCH, NPERM,
STRUCT_LO and STRUCT_HI are all environment overrides.

## Fixed after the first run attempt (20 November)

The first version hung with no output at all. The cause was mine: the file
search used `glob(..., recursive=True)` over five `**/*.csv` patterns across
the whole project tree, and it ran before any line was printed. On a shared
parallel filesystem holding read sets, assemblies and IQ-TREE run
directories, that stats every directory and takes many minutes, which is
indistinguishable from a hang.

Three changes:

- a banner and the parameters print before any filesystem work, so progress
  is visible immediately
- known locations are tried first (`09_Temporal_Signal/root_to_tip_LSDV_only.csv`
  and its siblings), and any fallback is a bounded `os.walk` at depth 3 that
  prunes the heavy directories
- the full column list of the chosen table is printed, so a column-detection
  failure can be diagnosed from the log rather than guessed at

Passing `TEMPORAL=` skips the search entirely and is the fastest route.
