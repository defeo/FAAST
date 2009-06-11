/*
	This file is part of the FAAST library.

	Copyright (c) 2009 Luca De Feo and Éric Schost.

	The most recent version of FAAST is available at http://www.lix.polytechnique.fr/~defeo/FAAST

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; see file COPYING. If not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/
/**
	\example using_infrastructure.c++
	This example illustrates how to use the FAAST::Field::switchContext(),
	FAAST::FieldElement::toInfrastructure(), FAAST::FieldPolynomial::toInfrastructure()
	and FAAST::Field::fromInfrastructure() methods.

	\warning This is for advanced use only, you shouldn't care about this on an
	ordinary use of the library.
*/

#include <faast.hpp>
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
