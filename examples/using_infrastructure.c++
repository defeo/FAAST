#include <Artin-Schreier.hpp>
using namespace FAAST;
void main() {
	Field<zz_p_algebra>& K = Field<zz_p_Algebra>::createField(2, 3);
	Field<zz_p_algebra>& L = K.ArtinSchreierExtension();

	/* Set the context to work in K. */
	K.switchContext();
	NTL::zz_pX a_polynomial, another_polynomial;

	/* Do something with a_polynomial and another_polynomial
	 * WARNING ! If you make library calls here, do not forget
	 * to set the context back with K.switchContext().
	 */
	...

	/* Convert a_polynomial to an element of K. */
	NTL::zz_pE an_element = NTL::to_zz_pE(a_polynomial);
	FieldElement<zz_p_Algebra> a = K.fromInfrastructure(an_element);

	/* Convert another_polynomial to an element of L. This makes sense since
	 * K and L have the same characteristic.
	 */
	L.switchContext();
	NTL::zz_pE another_element = NTL::to_zz_pE(another_polynomial);
	FieldElement<zz_p_Algebra> b = L.fromInfrastructure(another_element);
	
	/* Now work using the library. You don't have to worry about
	 * switching the context as the library does it for you.
	 */
	FieldElement<zz_p_Algebra> c = a >> L;
	c *= b;
	
	/* Extract the NTL representation out of c */
	NTL::zz_pE yet_another_element;
	c.toInfrastructure(yet_another_element);
	
	/* We call an \NTL function not available in the library
	 * and we convert the result back to a FieldElement
	 */
	L.switchContext();	// this is optional because toInfrastructure()
						// already switched the context.
	NTL:zz_pX yet_another_polynomial = NTL::rep(yet_another_element);
	NTL::LeftShift(yet_another_polynomial, yet_another_polynomial, 1);
	yet_another_element = NTL::to_zz_pE(yet_another_polynomial);
	c = L.fromInfrastructure(yet_another_element);
}
