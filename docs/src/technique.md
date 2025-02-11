# Techniques

## Homogenization


## Christoffel-Darboux Kernels


### Methods
```@docs
run_H1
run_H1CS
run_H2
run_H2CS
```

## Tips for modelling polynomial optimization problems
- When possible, explictly include a sphere/ball constraint (or multi-sphere/multi-ball constraints).
- When the feasible set is unbounded, try the homogenization technique introduced in [Homogenization for polynomial optimization with unbounded sets](https://link.springer.com/article/10.1007/s10107-022-01878-5).
- Scale the coefficients of the polynomial optimization problem to $[-1, 1]$.
- Scale the variables so that they take values in $[-1, 1]$ or $[0, 1]$.
- Try to include more (redundant) inequality constraints.