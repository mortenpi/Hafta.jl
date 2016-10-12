using Base.Test
import Hafta

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
