#include <Artin-Schreier.hpp>
using namespace AS;
void main() {
	Field<zz_p_algebra>& K = Field<zz_p_Algebra>::createField(2, 3);
	Field<zz_p_algebra>& L = K.ArtinSchreierExtension();

	/* Set the context to work in K */
	K.switchContext();
	zz_pX a_polynomial, another_polynomial;

	/* Do something with a_polynomial and another_polynomial
	 * WARNING ! If you make library calls here, do not forget
	 * to set the context back with K.switchContext()
	 */
	...

	/* Convert a_polynomial to an element of K */
	zz_pE an_element = to_zz_pE(a_polynomial);
	FieldElement<zz_p_Algebra> a = K.fromInfrastructure(an_element);

	/* Convert another_polynomial to an element of L. This makes sense since
	 * K and L have the same characteristic. */
	L.switchContext();
	zz_pE another_element = to_zz_pE(another_polynomial);
	FieldElement<zz_p_Algebra> b = L.fromInfrastructure(another_element);
}
