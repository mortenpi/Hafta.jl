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
