# Independent numeric algebra checks and regressions for quaternion equalities.
module QuaternionAlgebraAudit

using TSSOS, TSSOS.QuaternionSOS, DynamicPolynomials, Quaternions
using JuMP, Clarabel, LinearAlgebra, Random, Test

const QS = TSSOS.QuaternionSOS
const Iq = QuaternionF64(0, 1, 0, 0)

function evaluate_term(term, values, n)
    result = QuaternionF64(1)
    for (i, j) in zip(term[1], term[2])
        result *= values[i] * conj(values[j-n])
    end
    for i in term[3]
        result *= i <= n ? values[i] : conj(values[i-n])
    end
    return result
end

@testset "Independent quaternion algebra audit" begin
    @testset "Word reduction, adjoint, and cyclic real part" begin
        rng = MersenneTwister(517)
        for _ in 1:1000
            n = 3
            values = randn(rng, QuaternionF64, n)
            term = [UInt16[], UInt16[], rand(rng, UInt16(1):UInt16(2n), rand(rng, 0:9))]
            value = evaluate_term(term, values, n)
            @test evaluate_term(QS.standardterm(term, n), values, n) ≈ value atol=1e-9
            @test evaluate_term(QS.star(term, n), values, n) ≈ conj(value) atol=1e-9
            canonical = QS.canonical_qmono(term, n)
            @test real(evaluate_term(canonical, values, n)) ≈ real(value) atol=1e-9
            @test QS.canonical_qmono(canonical, n) == canonical
        end
    end

    @testset "Quaternion-to-real evaluation" begin
        rng = MersenneTwister(518)
        @ncpolyvar z[1:4]
        for _ in 1:50
            values = randn(rng, QuaternionF64, 2)
            word = rand(rng, 1:4, rand(rng, 1:6))
            coefficient = randn(rng, QuaternionF64)
            polynomial = coefficient * prod(z[i] for i in word)
            pop, x = quaternion_to_real([polynomial], z)
            coords = [real.(values); [imag_part(q)[j] for j in 1:3 for q in values]]
            expected = real(coefficient * prod(i <= 2 ? values[i] : conj(values[i-2]) for i in word))
            @test pop[1](x => coords) ≈ expected atol=1e-8
        end
    end

    @testset "Two-sided equality multiplier and adjoint" begin
        rng = MersenneTwister(519)
        for _ in 1:100
            values = randn(rng, QuaternionF64, 3)
            a, w, b = [QS.standardterm([UInt16[], UInt16[],
                       rand(rng, UInt16(1):UInt16(6), rand(rng, 0:3))], 3) for _ in 1:3]
            H, c = randn(rng, QuaternionF64, 2)
            components = QS.quaternion_right_product((real(H), imag_part(H)...), c)
            coefficient = QuaternionF64(components...)
            @test coefficient ≈ H*c atol=1e-10
            word = QS.qtermadd3(w, b, a, 3)
            represented = coefficient*evaluate_term(word, values, 3) +
                          conj(coefficient)*evaluate_term(QS.star(word, 3), values, 3)
            direct = conj(evaluate_term(a, values, 3))*H*c*
                     evaluate_term(w, values, 3)*evaluate_term(b, values, 3)
            @test real(represented) ≈ 2*real(direct) atol=1e-8
        end
    end

    @testset "Exact quaternion quadratic bounds" begin
        rng = MersenneTwister(513)
        @ncpolyvar z[1:4]
        for _ in 1:6
            A = randn(rng, QuaternionF64, 2, 2)
            Q = (A + A') / 2
            supp = [QS.standardterm([UInt16[], UInt16[], UInt16[j, i+2]], 2) for i in 1:2 for j in 1:2]
            coe = QuaternionF64[Q[i,j] for i in 1:2 for j in 1:2]
            a, d = real(Q[1,1]), real(Q[2,2])
            expected = min(0, (a+d-sqrt((a-d)^2+4abs2(Q[1,2]))) / 2)
            bound = qtssos([0, 1-z[1]*z[3]-z[2]*z[4]], z, 2, 1;
                           fsupp=supp, fcoe=coe, CS=false, QUIET=true, optimizer=Clarabel.Optimizer)
            @test bound ≈ expected atol=1e-6
        end
    end

    @testset "Quaternion equalities reject false bounds" begin
        @ncpolyvar q[1:2]
        for target in (QuaternionF64(1), Iq)
            # True optimum is exactly zero, attained at q = target.
            model = Model(Clarabel.Optimizer)
            bound = qtssos([0, q[1]-target], q, 1, 2;
                   model, numeq=1, qncnumeq=1, CS=false, QUIET=true)
            @test bound ≈ 0 atol=1e-6
            # Require an explicit finite false certificate, independently of
            # whether the original maximization reports unboundedness.
            fix(model[:lower], 1; force=true)
            @objective(model, Min, 0)
            optimize!(model)
            @test termination_status(model) == MOI.INFEASIBLE
        end
    end

    @testset "Quaternion multipliers enforce imaginary equalities" begin
        @ncpolyvar q[1:2]
        # q - conj(q) = 0 forces q to be real, so Re(i*q) is exactly zero.
        supp = [[UInt16[], UInt16[], UInt16[1]], [UInt16[], UInt16[], UInt16[2]]]
        coe = QuaternionF64[Iq/2, -Iq/2]
        for order in (1, 2), TS in (false, "block")
            bound = qtssos([0, 1-q[1]*q[2], q[1]-q[2]], q, 1, order;
                       fsupp=supp, fcoe=coe, numeq=1, qncnumeq=1,
                       CS=false, TS, QUIET=true, optimizer=Clarabel.Optimizer)
            @test bound ≈ 0 atol=1e-6
        end
    end

    @testset "Fixed quaternions match real coordinate constraints" begin
        @ncpolyvar q[1:2]
        @polyvar x[1:4]
        target = QuaternionF64(1, 2, -1, 0.5)
        coordinates = [real(target), imag_part(target)...]
        for (objective, real_objective, order, expected) in (
            (q[1]*q[2], sum(x.^2), 1, abs2(target)),
            (-q[1]*q[2], -sum(x.^2), 2, -abs2(target)),
            ((q[1]+q[2])/2, x[1], 1, real(target)),
        )
            real_bound = tssos([real_objective; x .- coordinates], x, order;
                               numeq=4, TS=false, QUIET=true, optimizer=Clarabel.Optimizer)[1]
            for TS in (false, "block")
                bound = qtssos([objective, q[1]-target], q, 1, order;
                               numeq=1, qncnumeq=1, CS=false, TS,
                               QUIET=true, optimizer=Clarabel.Optimizer)
                @test bound ≈ expected atol=1e-6
                @test bound ≈ real_bound atol=1e-6
            end
        end
    end
end

end
