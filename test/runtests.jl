if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

import Hafta

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

@testset "HFB" begin end
