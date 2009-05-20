/*
 * Obsolete (integrated in the main code)
 *
 * GF2hack.c++
 *
 *  Created on: Apr 29, 2009
 *      Author: defeo
 */

namespace AS {
	template <> FieldPolynomial<GF2_Algebra>&
	FieldPolynomial<GF2_Algebra>::operator=(const FieldPolynomial<GF2_Algebra>& e)
	throw() {
		// Hack for NTL
		if (parent_field && ! parent_field->isIsomorphic(*e.parent_field)) {
			repBase.kill();
			repExt.kill();
		}

		base = e.base;
		parent_field = e.parent_field;
		if (parent_field) {
			parent_field->switchContext();
			if (base) { repBase = e.repBase; repExt = 0; }
			else { repExt = e.repExt; repBase = 0; }
		} else {
			repBase.kill();
			repExt.kill();
		}
		return *this;
	}
}
