<!-- badges: start -->
[![R build status](https://github.com/jakubkala/QuiPTsim/workflows/R-CMD-check/badge.svg)](https://github.com/jakubkala/QuiPTsim/actions)
[![Codecov test coverage](https://codecov.io/gh/jakubkala/QuiPTsim/branch/master/graph/badge.svg)](https://codecov.io/gh/jakubkala/QuiPTsim?branch=master)
<!-- badges: end -->

# QuiPTsim
## Quick Permutation Test simulation

### TODO

- [X] add probabilities' vectors for both motifs and sequences
- [X] longer motifs: ds <= 6 ns <= 4
- [X] motifs' generation in each replication 
- [X] add function validate_motifs()
- [X] count multigrams for longer motifs
- [X] add sequences to attributes(dat)
- [X] tests

- [ ] add seqR counter

### Simulation details

* number of sequences: 3k positive / 3k negative
* replications: 50-100
* number of motifs: 1 - 3
* sequences' lengths: 10, 20, 40, 80
* alphabet: 20 elements: 4, 6, 8, 20
* various probability vectors

