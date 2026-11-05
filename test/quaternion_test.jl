module QuaternionTests

using Test, DynamicPolynomials, JuMP, Quaternions, Clarabel, TSSOS
using TSSOS.QuaternionSOS
const QS = TSSOS.QuaternionSOS
const opts = (; optimizer=Clarabel.Optimizer, QUIET=true)

@testset "Quaternion polynomial optimization" begin
    # Deliberately use local names other than q/n: the standalone code depended
    # on globals with these names in Main.
    @ncpolyvar z[1:4]
    norm1, norm2 = z[1]*z[3], z[2]*z[4]

    @testset "Algebra and conversion" begin
        @test quaternion_to_real([0, 1], z)[1] == [0, 1]
        pop, x = quaternion_to_real([-norm1, 1-norm1], z)
        @test pop[1] == -sum(x[[1,3,5,7]].^2)
        @test pop[2] == 1+pop[1]
        reversed_vars = reverse(z)
        supp, _ = QS.QPolys_info([z[1]*z[3]], reversed_vars, 2)
        @test supp[1] == [[UInt16[2], UInt16[4], UInt16[]]]
        term = [UInt16[1,2], UInt16[3,4], UInt16[]]
        @test QS.qreduce_unitnorm(term, 2; nb=1) == [UInt16[2], UInt16[4], UInt16[]]
        a = [UInt16[], UInt16[], UInt16[1]]
        b = [UInt16[], UInt16[], UInt16[2]]
        original = deepcopy(a)
        @test QS.qtermaddleft(a, b, 2)[3] == UInt16[3,2]
        QS.qtermadd(a, b, 2)
        @test a == original
        @test length(QS.get_qncbasis([1], 2, 2; conjubasis=true)) == 6
        @test QS.bfind === TSSOS.bfind
        @test QS.chordal_cliques! === TSSOS.chordal_cliques!
    end

    @testset "Dense, sparse, and unit-norm bounds" begin
        for ipart in (false, true), TS in (false, "block", "MD")
            bound = qtssos([-norm1, 1-norm1], z, 2, 1;
                           CS=false, TS, ipart, addrcons=false, opts...)
            @test bound ≈ -1 atol=1e-6
        end
        for CS in (false, "MF", "NC"), TS in (false, "block")
            @test qtssos([-norm1-norm2, 1-norm1, 1-norm2], z, 2, 1;
                         CS, TS, opts...) ≈ -2 atol=1e-6
        end
        @test qtssos([-norm1], z, 2, 1; nb=1, CS=false, opts...) ≈ -1 atol=1e-6
        @test qtssos([-norm1, 1-norm1], z, 2, 1;
                     numeq=1, CS=false, addrcons=false, opts...) ≈ -1 atol=1e-6
        @test qtssos([-norm1^2, 1-norm1], z, 2, 2;
                     CS="MF", addrcons=false, ipart=false, opts...) ≈ -1 atol=1e-6
        for ipart in (false, true)
            @test qtssos([norm1], z, 2, 1;
                         CS=false, TS="block", addrcons=false, ipart, opts...) ≈ 0 atol=1e-6
        end
        @test qtssos([-norm1, 1-norm1, z[2]], z, 2, 1;
                     CS=false, numeq=1, qncnumeq=1, opts...) ≈ -1 atol=1e-6
        @test qtssos([-norm1-norm2, 1-norm1, 1-norm2], z, 2, 1;
                     CS=false, addicons=true, opts...) ≈ -2 atol=1e-6
        @ncpolyvar v[1:6]
        @test qtssos([-v[1]*v[4], 1-v[1]*v[4], v[2]-v[3]], v, 3, 1;
                     numeq=1, qncnumeq=1, opts...) ≈ -1 atol=1e-6
        realpop, x = quaternion_to_real([-norm1, 1-norm1], z)
        realbound = tssos(realpop, x, 1; TS=false, opts...)[1]
        @test realbound ≈ -1 atol=1e-6
    end

    @testset "Explicit quaternion coefficients and model interface" begin
        @ncpolyvar w[1:2]
        @test qtssos([0], w, 1, 0; opts...) ≈ 0 atol=1e-6
        @test qtssos([0, 1-w[1]*w[2]], w, 1, 1; opts...) ≈ 0 atol=1e-6
        # Re(i*q) = -Im_i(q), whose minimum on the unit ball is -1.
        # Supply both adjoint terms of the Hermitian objective.
        fsupp = [[UInt16[], UInt16[], UInt16[1]], [UInt16[], UInt16[], UInt16[2]]]
        fcoe = QuaternionF64[QuaternionF64(0,0.5,0,0), QuaternionF64(0,-0.5,0,0)]
        bound = qtssos([0.5*(w[1]+w[2]), 1-w[1]*w[2]], w, 1, 1;
                       fsupp, fcoe, CS=false, opts...)
        @test bound ≈ -1 atol=1e-6
        model = Model(Clarabel.Optimizer)
        @test qtssos([-w[1]*w[2], 1-w[1]*w[2]], w, 1, 1;
                     model, QUIET=true) ≈ -1 atol=1e-6
        @test num_variables(model) > 0
        @test qtssos([-w[1]*w[2]], w, 1, 1; solve=false, QUIET=true) === nothing
        @test_throws ArgumentError qtssos([w[1]], w, 1, 1; fsupp, solve=false, QUIET=true)
        @test_throws ArgumentError qtssos([w[1]^4], w, 1, 1; solve=false, QUIET=true)
        @test_throws ArgumentError qtssos([w[1]], w[1:1], 1, 1; solve=false, QUIET=true)
        @test_throws ArgumentError qtssos([w[1]], w, 1, 1; model, opts...)
    end
end

end
