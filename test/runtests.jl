if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

import Hafta

@testset "Hafta tests" begin

@testset "HarmonicOscillator" begin
    @test Hafta.HarmonicOscillator.wmatrix_arridx_fromsorted(36,46,37,22) == 22907
    @test Hafta.HarmonicOscillator.wmatrix_arridx_fromsorted_recursive(36,46,37,22) == 22907

    # The empty_expression() function needs some tests, just to make sure that
    # it works properly with all versions.
    @test   Hafta.HarmonicOscillator.empty_expression( :() )
    @test ! Hafta.HarmonicOscillator.empty_expression( :(5,6) )
    @test ! Hafta.HarmonicOscillator.empty_expression( :([]) )
    # These fail at the moment because these expression are not Expr objects.
    # See the comment and TOOD at the function definition.
    #@test ! Hafta.HarmonicOscillator.empty_expression( :(x) )
    #@test ! Hafta.HarmonicOscillator.empty_expression( :(5) )

    # TODO: the basis function definitions should be switched to types probably.
    # Then we'd be storing the polynomial as well and we could check parity with
    # that.
    for i in 0:2:100 @test Hafta.HarmonicOscillator.Psi_parity(i) === 1 end
    for i in 1:2:101 @test Hafta.HarmonicOscillator.Psi_parity(i) === -1 end
end

@testset "HarmonicSystems" begin
    let s = Hafta.Harmonic1DFermionSystem(3, -1.0)
        @test length(s) === 6
        @test s.vcoef === -1.0
        for i in 1:length(s)
            @test Hafta.H0(s, i, i) === div(i-1, 2) + 0.5
        end
        @test Hafta.HarmonicSystems.nshells(s) === 3
    end
end

@testset "HartreeFock" begin end

@testset "HFB" begin
    # The HFB module uses more optimized routines now to calculate the Vbar and
    # Γ and Δ fields. This compares them to the older and slower routines.
    Vbar_simple(s::Hafta.HFB.HFBSystem, i,j,k,l) = Hafta.V(s.system, i,j,k,l) - Hafta.V(s.system, i,j,l,k)

    s = Hafta.HFB.HFBSystem(Hafta.Harmonic2DSystem(3, -1.0))
    @test Hafta.HFB.Vbar(s, 1,1,1,1) == Vbar_simple(s, 1,1,1,1)
    @test Hafta.HFB.Vbar(s, 1,2,1,2) == Vbar_simple(s, 1,2,1,2)
    @test Hafta.HFB.Vbar(s, 1,2,2,1) == Vbar_simple(s, 1,2,2,1)

    function gamma_delta_slow(system::Hafta.HFB.HFBSystem, rho::Matrix, kappa::Matrix)
        M = length(system)
        gamma = zeros(Float64, (M,M))
        delta = zeros(Float64, (M,M))
        for i=1:M,j=1:M,k=1:M,l=1:M
            gamma[i,j] += rho[k,l] * Vbar_simple(system, i,k,j,l)
            delta[i,j] += 0.5*kappa[k,l] * Vbar_simple(system, i,j,k,l)
        end
        gamma,delta
    end

    function are_ΓΔ_equal(ΓΔ1, ΓΔ2)
        Γ1, Δ1 = ΓΔ1
        Γ2, Δ2 = ΓΔ2
        Γ1 ≈ Γ2 && Δ1 ≈ Δ2
    end

    let M = length(s), rho = rand(Float64, (M,M)), kappa = rand(Float64, (M,M))
        ΓΔ_slow = @time gamma_delta_slow(s, rho, kappa)
        ΓΔ = @time Hafta.HFB.gamma_delta(s, rho, kappa)
        @test are_ΓΔ_equal(ΓΔ_slow, ΓΔ)
    end
end

include("qrpa.jl")

end
