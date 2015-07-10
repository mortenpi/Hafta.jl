using Base.Test
import Hafta

@test Hafta.HarmonicOscillator.wmatrix_arridx_fromsorted(36,46,37,22) == 22907
@test Hafta.HarmonicOscillator.wmatrix_arridx_fromsorted_recursive(36,46,37,22) == 22907
