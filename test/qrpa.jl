if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

import Hafta

@testset "QRPA" begin
    hfbi, op, solutions = Hafta.QRPA.get_qrpa_values(Hafta.HarmonicSystems.Harmonic2DSystem, 3, 2, -6.1; nev=5, verbose=true)
    @test length(solutions) === 5
    @test last(hfbi.es) ≈ 0.7445101513069428
    for s in solutions
        @test imag(s.energy) ≈ 0.0
    end
    local energies = [real(s.energy) for s in solutions]
    @test maximum(energies) ≈ 1.3037555525299303
end
