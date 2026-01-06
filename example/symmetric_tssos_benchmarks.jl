include("symmetric_tssos_benchmark_tools.jl");


#####################################################
## Example 1D ring Ising quartic on Laplacian ball ##
#####################################################

n=6;
X, f, g, G, r = build_ising_quartic(n; a=1.0, b=0.5, c=-1.0, d=0.2);
pop=[f,g];


# dense
opt0,data0 = tssos(pop, X, r; numeq = 0, TS=false);

# sign symmetry
opt1,sol1,data1 = tssos(pop, X, r; numeq = 0, TS="signsymmetry");
opt_higher1,sol1,data1 = tssos(data1; TS="signsymmetry");

# tssos (max chordal extension)
opt2,sol2,data2 = tssos(pop, X, r; numeq = 0, TS="block", QUIET=false);
opt_higher2,sol2,data2 = tssos(data2; TS="block");

# symmetry
opt3,data3 = tssos_symmetry(pop, X, r, G; numeq=0, SymmetricConstraint=true);

# symmetry + tssos (max chordal extension, WITH diagonal squares)
opt4,data4 = tssos_symmetry(pop, X, r, G; numeq=0, DiagSquare=true, SymmetricConstraint=true, TS="block");
opt_higher4,data4 = tssos_symmetry(data4; TS="block");

# symmetry + tssos (max chordal extension, WITHOUT diagonal squares)
opt5,data5 = tssos_symmetry(pop, X, r, G; numeq=0, DiagSquare=false, SymmetricConstraint=true, TS="block");
opt_higher5,data5 = tssos_symmetry(data5; TS="block");



####################################################
## Example 2D torus grid quartic on discrete cube ##
####################################################

p,q=2,3;
n=p*q;
X, f, g, G, r = build_grid_quartic(p, q; a=1.0, b=0.5, c=-1.0, d=0.2, opt=:spin);
pop=[f,g];

# opt (:smooth, :magnet, :spin)
# for opt=:spin, treat g as an equality and set numeq=1


# dense
opt0,data0 = tssos(pop, X, r; numeq = 1, TS=false);

# sign symmetry
opt1,sol1,data1 = tssos(pop, X, r; numeq = 1, TS="signsymmetry");
opt_higher1,sol1,data1 = tssos(data1; TS="signsymmetry");

# tssos (max chordal extension)
opt2,sol2,data2 = tssos(pop, X, r; numeq = 1, TS="block");
opt_higher2,sol2,data2 = tssos(data2; TS="block");

# symmetry
opt3,data3 = tssos_symmetry(pop, X, r, G; numeq=1, SymmetricConstraint=true);

# symmetry + tssos (max chordal extension, WITH diagonal squares)
opt4,data4 = tssos_symmetry(pop, X, r, G; numeq=1, DiagSquare=true, SymmetricConstraint=true, TS="block");
opt_higher4,data4 = tssos_symmetry(data4; TS="block");

# symmetry + tssos (max chordal extension, WITHOUT diagonal squares)
opt5,data5 = tssos_symmetry(pop, X, r, G; numeq=1, DiagSquare=false, SymmetricConstraint=true, TS="block");
opt_higher5,data5 = tssos_symmetry(data5; TS="block");



#######################################
## Comparison with Symmetric Chordal ##
#######################################

n=6;
X, f, G, r = build_symmetric_quartic(n);


# tssos (max chordal ext)
opt,sol,data = tssos([f], X, d; numeq=0, TS="block");
opt_higher,data = tssos(data; TS="block");

# tssos (approx chordal ext)
opt,sol,data = tssos([f], X, d; numeq=0, TS="MD");
opt_higher,data_higher = tssos(data; TS="MD");

# symmetry + tssos (max chordal ext, WITH diagonal squares)
opt,data = tssos_symmetry([f], X, d, G; numeq=0, TS="block", DiagSquare=true);
opt_higher,data = tssos_symmetry(data; TS="block");

# symmetry + tssos (approx chordal ext, WITH diagonal squares)
opt,data = tssos_symmetry([f], X, d, G; numeq=0, TS="MD", DiagSquare=true);
opt_higher,data_higher = tssos_symmetry(data; TS="MD");

# symmetry + tssos (max chordal ext, WITHOUT diagonal squares)
opt,data = tssos_symmetry([f], X, d, G; numeq=0, TS="block", DiagSquare=false);
opt_higher,data = tssos_symmetry(data; TS="block");

# symmetry + tssos (approx chordal ext, WITHOUT diagonal squares)
opt,data = tssos_symmetry([f], X, d, G; numeq=0, TS="MD", DiagSquare=false);
opt_higher,data = tssos_symmetry(data; TS="MD");