# TSSOS

[TSSOS](https://github.com/wangjie212/TSSOS) is a sparse polynomial optimization package based on the sparsity adapted moment-SOS hierarchies, which can fully exploit the sparsity in the problem data including correlative (variable) sparsity and term sparsity.

---

### Authors

- [Jie Wang](https://wangjie212.github.io/jiewang), Academy of Mathematics and Systems Science, Chinese Academy of Sciences.

### Installation

[TSSOS](https://github.com/wangjie212/TSSOS) is simply installed by running

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

1. [TSSOS: A Moment-SOS hierarchy that exploits term sparsity](https://arxiv.org/abs/1912.08899), Jie Wang, Victor Magron and Jean B. Lasserre, 2020.
2. [Chordal-TSSOS: a moment-SOS hierarchy that exploits term sparsity with chordal extension](https://arxiv.org/abs/2003.03210), Jie Wang, Victor Magron and Jean B. Lasserre, 2020.
3. [CS-TSSOS: Correlative and term sparsity for large-scale polynomial optimization](https://arXiv:2005.02828), Jie Wang, Victor Magron, Jean B. Lasserre and Ngoc H. A. Mai, 2020.
4. [Exploiting term sparsity in Noncommutative Polynomial Optimization](https://arxiv.org/abs/2010.06956), Jie Wang and Victor Magron, 2020.
