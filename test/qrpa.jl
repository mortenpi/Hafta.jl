if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

import Hafta

@testset "QRPA" begin
    hfbi, op, solutions = Hafta.QRPA.get_qrpa_values(Hafta.HarmonicSystems.Harmonic2DSystem, 3, 2, -6.1; nev=10, verbose=true)
    @test length(solutions) === 10
    @test last(hfbi.es) ≈ 0.7445101513069428
    for s in solutions
        @test imag(s.energy) ≈ 0.0
    end
    local energies = [abs(real(s.energy)) for s in solutions]
    @test minimum(energies) ≈ 0.0031147182602445644
end
