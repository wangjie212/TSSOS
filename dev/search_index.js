var documenterSearchIndex = {"docs":
[{"location":"ncpop/#Noncommutative-polynomial-optimization","page":"Noncommutative Polynomial Optimization","title":"Noncommutative polynomial optimization","text":"","category":"section"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"TSSOS supports noncommutative polynomial (eigenvalue or trace) optimization involving noncommuting variables. See also NCTSSOS.","category":"page"},{"location":"ncpop/","page":"Noncommutative Polynomial Optimization","title":"Noncommutative Polynomial Optimization","text":"nctssos_first\nnctssos_higher!\ncs_nctssos_first\ncs_nctssos_higher!","category":"page"},{"location":"ncpop/#TSSOS.NCTSSOS.nctssos_first","page":"Noncommutative Polynomial Optimization","title":"TSSOS.NCTSSOS.nctssos_first","text":"opt,data = nctssos_first(f::Polynomial{false, T} where T<:Number, x::Vector{PolyVar{false}};\n    newton=true, reducebasis=true, TS=\"block\", obj=\"eigen\", merge=false, md=3, solve=true, QUIET=false)\n\nCompute the first step of the NCTSSOS hierarchy for unconstrained noncommutative polynomial optimization. If newton=true, then compute a monomial basis by the Newton chip method. If reducebasis=true, then remove monomials from the monomial basis by diagonal inconsistency. If TS=\"block\", use maximal chordal extensions; if TS=\"MD\", use approximately smallest chordal extensions. If obj=\"eigen\", minimize the eigenvalue; if obj=\"trace\", then minimize the trace. If merge=true, perform the PSD block merging. Return the optimum and other auxiliary data.\n\nArguments\n\nf: the objective function for unconstrained noncommutative polynomial optimization.\nx: the set of noncommuting variables.\nmd: the tunable parameter for merging blocks.\n\n\n\n\n\n","category":"function"},{"location":"ncpop/#TSSOS.NCTSSOS.nctssos_higher!","page":"Noncommutative Polynomial Optimization","title":"TSSOS.NCTSSOS.nctssos_higher!","text":"opt,data = nctssos_higher!(data, TS=\"block\", merge=false, md=3, QUIET=false)\n\nCompute higher steps of the NCTSSOS hierarchy. Return the optimum and other auxiliary data.\n\n\n\n\n\n","category":"function"},{"location":"ncpop/#TSSOS.NCTSSOS.cs_nctssos_first","page":"Noncommutative Polynomial Optimization","title":"TSSOS.NCTSSOS.cs_nctssos_first","text":"opt,data = cs_nctssos_first(pop, x, d; numeq=0, CS=\"MF\", TS=\"block\", merge=false, md=3,\nQUIET=false, obj=\"eigen\", solve=true)\n\nCompute the first step of the CS-NCTSSOS hierarchy for constrained noncommutative polynomial optimization with relaxation order d. Return the optimum and other auxiliary data.\n\nArguments\n\npop: the vector of the objective function, inequality constraints, and equality constraints.\nx: the set of noncommuting variables.\nd: the relaxation order of the moment-SOHS hierarchy.\nnumeq: the number of equality constraints.\n\n\n\n\n\nopt,data = cs_nctssos_first(supp::Vector{Vector{Vector{UInt16}}}, coe, n::Int, d::Int; numeq=0,\nCS=\"MF\", TS=\"block\", merge=false, md=3, QUIET=false, obj=\"eigen\", solve=true)\n\nCompute the first step of the CS-NCTSSOS hierarchy for constrained noncommutative polynomial optimization with relaxation order d. Here the polynomial optimization problem is defined by supp and coe, corresponding to the supports and coeffients of pop respectively. Return the optimum and other auxiliary data.\n\nArguments\n\nsupp: the supports of the polynomial optimization problem.\ncoe: the coeffients of the polynomial optimization problem.\nd: the relaxation order of the moment-SOHS hierarchy.\nnumeq: the number of equality constraints.\n\n\n\n\n\n","category":"function"},{"location":"ncpop/#TSSOS.NCTSSOS.cs_nctssos_higher!","page":"Noncommutative Polynomial Optimization","title":"TSSOS.NCTSSOS.cs_nctssos_higher!","text":"opt,data = cs_nctssos_higher!(data; TS=\"block\", QUIET=false, merge=false, md=3, solve=true)\n\nCompute higher steps of the CS-NCTSSOS hierarchy. Return the optimum and other auxiliary data.\n\n\n\n\n\n","category":"function"},{"location":"#TSSOS","page":"Home","title":"TSSOS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TSSOS is a sparse polynomial optimization package based on the sparsity adapted moment-SOS hierarchies, which can fully exploit the sparsity in the problem data including correlative (variable) sparsity and term sparsity.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Jie Wang, Academy of Mathematics and Systems Science, Chinese Academy of Sciences.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TSSOS is simply installed by running","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/wangjie212/TSSOS","category":"page"},{"location":"#Related-packages","page":"Home","title":"Related packages","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"DynamicPolynomials: Polynomial definition\nMultivariatePolynomials: Polynomial manipulation\nNCTSSOS: Noncommutative polynomial optimization\nChordalGraph: Chordal graphs and chordal extentions\nSparseJSR: Computing joint spetral radius","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TSSOS: A Moment-SOS hierarchy that exploits term sparsity, Jie Wang, Victor Magron and Jean B. Lasserre, 2020.\nChordal-TSSOS: a moment-SOS hierarchy that exploits term sparsity with chordal extension, Jie Wang, Victor Magron and Jean B. Lasserre, 2020.\nCS-TSSOS: Correlative and term sparsity for large-scale polynomial optimization, Jie Wang, Victor Magron, Jean B. Lasserre and Ngoc H. A. Mai, 2020.\nExploiting term sparsity in Noncommutative Polynomial Optimization, Jie Wang and Victor Magron, 2020.","category":"page"},{"location":"spop/#Polynomial-optimization","page":"Polynomial Optimization","title":"Polynomial optimization","text":"","category":"section"},{"location":"spop/","page":"Polynomial Optimization","title":"Polynomial Optimization","text":"tssos_first\ntssos_higher!\ncs_tssos_first\ncs_tssos_higher!\nrefine_sol","category":"page"},{"location":"spop/#TSSOS.tssos_first","page":"Polynomial Optimization","title":"TSSOS.tssos_first","text":"opt,sol,data = tssos_first(f, x; nb=0, newton=true, reducebasis=false, TS=\"block\", merge=false,\nmd=3, feasible=false, solver=\"Mosek\", QUIET=false, solve=true, MomentOne=false, Gram=false, solution=false, tol=1e-4)\n\nCompute the first TS step of the TSSOS hierarchy for unconstrained polynomial optimization. If newton=true, then compute a monomial basis by the Newton polytope method. If reducebasis=true, then remove monomials from the monomial basis by diagonal inconsistency. If TS=\"block\", use maximal chordal extensions; if TS=\"MD\", use approximately smallest chordal extensions.  If merge=true, perform the PSD block merging.  If feasible=true, then solve the feasibility problem. If solve=false, then do not solve the SDP. If Gram=true, then output the Gram matrix. If MomentOne=true, add an extra first order moment matrix to the moment relaxation.\n\nInput arguments\n\nf: objective\nx: POP variables\nnb: number of binary variables in x\nmd: tunable parameter for merging blocks\nQUIET: run in the quiet mode (true, false)\ntol: relative tolerance to certify global optimality\n\nOutput arguments\n\nopt: optimum\nsol: (near) optimal solution (if solution=true)\ndata: other auxiliary data \n\n\n\n\n\nopt,sol,data = tssos_first(pop, x, d; nb=0, numeq=0, quotient=true, basis=[],\nreducebasis=false, TS=\"block\", merge=false, md=3, solver=\"Mosek\", QUIET=false, solve=true,\nMomentOne=false, Gram=false, solution=false, tol=1e-4)\n\nCompute the first TS step of the TSSOS hierarchy for constrained polynomial optimization. If quotient=true, then exploit the quotient ring structure defined by the equality constraints. If merge=true, perform the PSD block merging.  If solve=false, then do not solve the SDP. If Gram=true, then output the Gram matrix. If MomentOne=true, add an extra first order moment matrix to the moment relaxation.\n\nInput arguments\n\npop: vector of the objective, inequality constraints, and equality constraints\nx: POP variables\nd: relaxation order\nnb: number of binary variables in x\nnumeq: number of equality constraints\nTS: type of term sparsity (\"block\", \"MD\", \"MF\", false)\nmd: tunable parameter for merging blocks\nnormality: impose the normality condtions (true, false)\nQUIET: run in the quiet mode (true, false)\ntol: relative tolerance to certify global optimality\n\nOutput arguments\n\nopt: optimum\nsol: (near) optimal solution (if solution=true)\ndata: other auxiliary data \n\n\n\n\n\n","category":"function"},{"location":"spop/#TSSOS.tssos_higher!","page":"Polynomial Optimization","title":"TSSOS.tssos_higher!","text":"opt,sol,data = tssos_higher!(data; TS=\"block\", merge=false, md=3, QUIET=false, solve=true,\nMomentOne=false, solution=false, tol=1e-4)\n\nCompute higher TS steps of the TSSOS hierarchy.\n\n\n\n\n\n","category":"function"},{"location":"spop/#TSSOS.cs_tssos_first","page":"Polynomial Optimization","title":"TSSOS.cs_tssos_first","text":"opt,sol,data = cs_tssos_first(pop, x, d; nb=0, numeq=0, CS=\"MF\", cliques=[], TS=\"block\", merge=false, md=3, solver=\"Mosek\", QUIET=false, solve=true, solution=false,\nGram=false, MomentOne=false, Mommat=false, tol=1e-4)\n\nCompute the first TS step of the CS-TSSOS hierarchy for constrained polynomial optimization. If merge=true, perform the PSD block merging.  If solve=false, then do not solve the SDP. If Gram=true, then output the Gram matrix. If Mommat=true, then output the moment matrix. If MomentOne=true, add an extra first order moment matrix to the moment relaxation.\n\nInput arguments\n\npop: vector of the objective, inequality constraints, and equality constraints\nx: POP variables\nd: relaxation order\nnb: number of binary variables in x\nnumeq: number of equality constraints\nCS: method of chordal extension for correlative sparsity (\"MF\", \"MD\", \"NC\", false)\ncliques: the set of cliques used in correlative sparsity\nTS: type of term sparsity (\"block\", \"MD\", \"MF\", false)\nmd: tunable parameter for merging blocks\nnormality: impose the normality condtions (true, false)\nQUIET: run in the quiet mode (true, false)\ntol: relative tolerance to certify global optimality\n\nOutput arguments\n\nopt: optimum\nsol: (near) optimal solution (if solution=true)\ndata: other auxiliary data \n\n\n\n\n\nopt,sol,data = cs_tssos_first(supp::Vector{Vector{Vector{UInt16}}}, coe, n, d; nb=0, numeq=0, CS=\"MF\", cliques=[], TS=\"block\", \nmerge=false, md=3, QUIET=false, solver=\"Mosek\", solve=true, solution=false, Gram=false, MomentOne=false, Mommat=false, tol=1e-4)\n\nCompute the first TS step of the CS-TSSOS hierarchy for constrained polynomial optimization.  Here the polynomial optimization problem is defined by supp and coe, corresponding to the supports and coeffients of pop respectively.\n\n\n\n\n\nopt,sol,data = cs_tssos_first(pop, z, n, d; nb=0, numeq=0, CS=\"MF\", cliques=[], TS=\"block\", ipart=true, merge=false, md=3, \nsolver=\"Mosek\", QUIET=false, solve=true, Gram=false, Mommat=false, MomentOne=false)\n\nCompute the first TS step of the CS-TSSOS hierarchy for constrained complex polynomial optimization.  If merge=true, perform the PSD block merging.  If ipart=false, then use the real moment-HSOS hierarchy. If solve=false, then do not solve the SDP. If Gram=true, then output the Gram matrix. If Mommat=true, then output the moment matrix. If MomentOne=true, add an extra first order moment matrix to the moment relaxation.\n\nInput arguments\n\npop: vector of the objective, inequality constraints, and equality constraints\nz: CPOP variables and their conjugate\nn: number of CPOP variables\nd: relaxation order\nnb: number of unit-norm variables in x\nnumeq: number of equality constraints\nCS: method of chordal extension for correlative sparsity (\"MF\", \"MD\", false)\ncliques: the set of cliques used in correlative sparsity\nTS: type of term sparsity (\"block\", \"MD\", \"MF\", false)\nmd: tunable parameter for merging blocks\nnormality: normal order\nQUIET: run in the quiet mode (true, false)\n\nOutput arguments\n\nopt: optimum\ndata: other auxiliary data \n\n\n\n\n\nopt,sol,data = cs_tssos_first(supp::Vector{Vector{Vector{Vector{UInt16}}}}, coe::Vector{Vector{ComplexF64}},\nn, d; nb=0, numeq=0, CS=\"MF\", cliques=[], TS=\"block\", ipart=true, merge=false, md=3, solver=\"Mosek\",\nQUIET=false, solve=true, Gram=false, Mommat=false, MomentOne=false)\n\nCompute the first TS step of the CS-TSSOS hierarchy for constrained complex polynomial optimization.  Here the complex polynomial optimization problem is defined by supp and coe, corresponding to the supports and coeffients of pop respectively.\n\n\n\n\n\n","category":"function"},{"location":"spop/#TSSOS.cs_tssos_higher!","page":"Polynomial Optimization","title":"TSSOS.cs_tssos_higher!","text":"opt,sol,data = cs_tssos_higher!(data; TS=\"block\", merge=false, md=3, QUIET=false, solve=true,\nsolution=false, Gram=false, Mommat=false, MomentOne=false)\n\nCompute higher TS steps of the CS-TSSOS hierarchy.\n\n\n\n\n\nopt,sol,data = cs_tssos_higher!(data; TS=\"block\", merge=false, md=3, QUIET=false, solve=true,\nsolution=false, Gram=false, Mommat=false, MomentOne=false)\n\nCompute higher TS steps of the CS-TSSOS hierarchy.\n\n\n\n\n\n","category":"function"},{"location":"spop/#TSSOS.refine_sol","page":"Polynomial Optimization","title":"TSSOS.refine_sol","text":"ref_sol,flag = refine_sol(opt, sol, data, QUIET=false, tol=1e-4)\n\nRefine the obtained solution by a local solver. Return the refined solution, and flag=0 if global optimality is certified, flag=1 otherwise.\n\n\n\n\n\n","category":"function"}]
}
