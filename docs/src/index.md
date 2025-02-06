# TSSOS

[TSSOS](https://github.com/wangjie212/TSSOS) aims to provide a user-friendly and efficient tool for solving optimization problems with polynomials, which is based on the structured moment-SOS hierarchy.

### Authors

- [Jie Wang](https://wangjie212.github.io/jiewang), Academy of Mathematics and Systems Science, Chinese Academy of Sciences.

- [Victor Magron](https://homepages.laas.fr/vmagron), Laboratoire d'Architecture et Analyse des Systèmes, CNRS.

### Installation

[TSSOS](https://github.com/wangjie212/TSSOS) could be installed by running

```julia
pkg> add https://github.com/wangjie212/TSSOS
```

### Related packages

- [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl): Polynomial definition
- [MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl): Polynomial manipulation
- [NCTSSOS](https://github.com/wangjie212/NCTSSOS): Noncommutative polynomial optimization
- [ChordalGraph](https://github.com/wangjie212/ChordalGraph): Chordal graphs and chordal extentions
- [SparseJSR](https://github.com/wangjie212/SparseJSR): Computing joint spetral radius

### References

1. [TSSOS: A Moment-SOS hierarchy that exploits term sparsity](https://epubs.siam.org/doi/abs/10.1137/19M1307871), Jie Wang, Victor Magron, and Jean B. Lasserre, 2021.
2. [Chordal-TSSOS: a moment-SOS hierarchy that exploits term sparsity with chordal extension](https://epubs.siam.org/doi/10.1137/20M1323564), Jie Wang, Victor Magron, and Jean B. Lasserre, 2021.
3. [CS-TSSOS: Correlative and term sparsity for large-scale polynomial optimization](https://dl.acm.org/doi/abs/10.1145/3569709), Jie Wang, Victor Magron, Jean B. Lasserre, and Ngoc H. A. Mai, 2022.
4. [TSSOS: a Julia library to exploit sparsity for large-scale polynomial optimization](https://arxiv.org/abs/2103.00915), Victor Magron and Jie Wang, 2021.
