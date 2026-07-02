# CRAN comments for oalasso 1.0.0

## Test environments

* local: macOS (Apple Silicon), R 4.5.0
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1),
  windows-latest (R release), macos-latest (R release)
* win-builder: R Under development (unstable) (2026-06-29 r90199 ucrt)

## R CMD check results

Local and GitHub Actions: 0 errors | 0 warnings | 0 notes.
win-builder (R-devel): 0 errors | 0 warnings | 1 note —
"checking CRAN incoming feasibility": New submission; possibly misspelled
words in DESCRIPTION (Balde, Ertefaie, Lefebvre, Shortreed).

## Comments

* This is a new submission.
* The flagged words are the surnames of the authors of the implemented
  methods (Shortreed & Ertefaie 2017 <doi:10.1111/biom.12679>; Balde,
  Yang & Lefebvre 2023 <doi:10.1111/biom.13683>).
* The suggested package 'psAve' is on CRAN (submitted and accepted prior
  to this submission); all uses are conditional via requireNamespace().
* All other Suggests ('MatchIt', 'WeightIt', 'rpart') are used
  conditionally; core functionality runs with Imports ('glmnet',
  'cobalt') only.
