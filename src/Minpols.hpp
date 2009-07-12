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
namespace FAAST {
/****************** Minimal polynomials ******************/
	/* All the minimal polynomials up to the field F.
	 *		res = minimalPolynomials(F);
	 * res[0] contains minimalPolynomials(F), res[1] contains
	 * minimalPolynomials(F.overfield) and so on.
	 */
	template <class T> void FieldElement<T>::minimalPolynomials(
	const Field<T>& F, vector<FieldPolynomial<T> >& res)
	const throw(NotASubFieldException, NotSupportedException) {
		const Field<T>* G = parent_field->stem;
		if (!F.isSubFieldOf(*G)) throw NotASubFieldException();
		if (F.isPrimeField() && !F.isBaseField())
			throw NotSupportedException();

		long levels = G->height - F.height + 1;
		res.resize(levels);
		// minimal polynomial is  (X-this)
		levels--;
		res[levels] = -(*this);
		res[levels].setCoeff(1);
		res[levels] >>= *G;
		// first loop
		// go down as long as this is a member of the subfield
		try {
			while (G != F.stem) {
				res[levels - 1] =
					res[levels] >> *(G->subfield);
				levels--;
				G = G->subfield;
			}
		} catch (IllegalCoercionException e) {}
		// second loop
		// go down by galois conjugation
		while (G != F.stem) {
			FieldPolynomial<T> tmp = res[levels];
			res[levels - 1] = res[levels];
			for (BigInt i = 1 ; i < G->p ; i++) {
				tmp.self_frobenius(G->subfield->d);
				res[levels -1] *= tmp;
			}
			res[levels - 1] >>= *(G->subfield);
			G = G->subfield;
			levels--;
		}
		// move out of the stem
		res[0] >>= F;
	}

	/* The a-affine minimal polynomial over the field F,
	 * That is the minimum degree polynomial P of F[X] such that
	 * 		P(this) = a.
	 * The optional parameter minpols must contain either the result of
	 * this.minimalPolynomials(F,v), or must be an empty vector, in which
	 * case it is filled with the result of minimalPolynomials(F,v).
	 *
	 * throws: NotASubFieldException if F is not a subfield of this.parent
	 * throws: NoSuchPolynomialException if a is not in the field
	 * 			generated by this
	 * throws: NotSupportedException if F is not part of an Artin-Schreier tower
	 */
	template <class T> FieldPolynomial<T>
	FieldElement<T>::affineMinimalPolynomial(const Field<T>& F,
	const FieldElement<T>& a, vector<FieldPolynomial<T> >& minpols)
	const throw(NotASubFieldException, NoSuchPolynomialException,
	NotSupportedException, BadParametersException) {
		const Field<T>* G = parent_field->stem;
		const Field<T>* H = a.parent_field->stem;
		if (!F.isSubFieldOf(*G)) throw NotASubFieldException();

		if (F.isPrimeField() && !F.isBaseField())
			throw NotSupportedException();
		if (!H->isSubFieldOf(*G)) throw NoSuchPolynomialException();

		// if a is in F, simply return the constant polynomial
		if (G->isSubFieldOf(F)) return a >> F;

		// first loop :
		// go down as long as this is a member of the subfield
		FieldPolynomial<T> tmpa = a;
		tmpa >>= *H;
		FieldElement<T> tmpx = *this >> *G;
		try {
			while (G != F.stem) {
				tmpx >>= *(G->subfield);
				try {
					if (H == G) {
						tmpa >>= *(H->subfield);
						H = H->subfield;
					}
				} catch (IllegalCoercionException e) {
					throw NoSuchPolynomialException();
				}
				G = G->subfield;
			}
		} catch (IllegalCoercionException e) {}
		// get the list of minimal polynomials of this
		if (minpols.size() == 0)
			tmpx.minimalPolynomials(F, minpols);
		else if (long(minpols.size()) < G->height - F.height + 1)
			throw BadParametersException("minpols too small.");

		// second loop :
		// go down as long as a is a member of the subfield
		try {
			while (H != F.stem) {
				tmpa >>= *(H->subfield);
				H = H->subfield;
			}
			// if we arrived at the bottom, the interpolating
			// polynomial is simply a
			return FieldPolynomial<T>(tmpa) >> F;
		} catch (IllegalCoercionException e) {}

		// The interpolation begins, but we act like if a = 1
		// till G equals H
		FieldPolynomial<T> M = minpols[0].derivative();
		FieldPolynomial<T> res = tmpx.evaluate(M, minpols).inv();

		// third loop : interpolation
		long i = minpols.size() - 1;
		while (G != F.stem) {
			// at this point we take a into account
			if (G == H) res *= tmpa;
			res *= (minpols[i-1] >> *G) / minpols[i];
			FieldPolynomial<T> tmp = res;
			for (BigInt j = 1 ; j < G->p ; j++) {
				tmp.self_frobenius(G->subfield->d);
				res += tmp;
			}
			res >>= *(G->subfield);
			G = G->subfield;
			i--;
		}

		return res >> F;
	}

	/* Returns P(this).
	 *
	 * The optional parameter minpols must contain either the result of
	 * this.minimalPolynomials(P.parent,v), or must be an empty vector,
	 * in which case it is filled with the result of
	 * minimalPolynomials(P.parent,v).
	 *
	 * throws : IllegalCoercionException if this cannot be coerced to P.parent
	 * 			nor can P be coerced to this.parent
	 */
	template <class T> FieldElement<T>
	FieldElement<T>::evaluate(const FieldPolynomial<T>& P,
	vector<FieldPolynomial<T> >& minpols)
	const throw(IllegalCoercionException, BadParametersException) {
		const Field<T>* G = parent_field->stem;
		const Field<T>* F = P.parent_field;
		// first case : this must go up
		if (F->isOverFieldOf(*G)) {
			FieldElement<T> x = *this >> *F;
			// Horner rule
			FieldElement<T> coeff, res = F->zero();
			if (P != 0) P.getCoeff(P.degree(), res);
			for (long i = P.degree()-1 ; i >= 0 ; i--) {
				res *= x;
				P.getCoeff(i, coeff);
				res += coeff;
			}
			return res;
		}
		// second case : this must go down
		else if (F->isSubFieldOf(*G)) {
			// first loop :
			// go down as long as this is a member of the subfield
			FieldElement<T> tmp = *this;
			try {
				while (*G != *(F->stem)) {
					tmp >>= *(G->subfield);
					G = G->subfield;
				}
			} catch (IllegalCoercionException e) {}

			// get the minimal polynomials if needed
			if (minpols.size() == 0)
				tmp.minimalPolynomials(*F, minpols);
			else if (long(minpols.size()) < G->height - F->height + 1)
				throw BadParametersException("minpols too small.");

			// second loop :
			// go up reducing by the minimal polynomials
			FieldPolynomial<T> reduced = P;
			long i = 0;
			try {
				F = F->stem;
				while (*F != *G) {
					reduced %= minpols[i];
					reduced >>= *(F->overfield);
					F = F->overfield;
					i++;
				}
			} catch (IllegalCoercionException e) {
				throw BadParametersException("Bad polynomials in minpols.");
			}
			// Horner rule
			FieldElement<T> res = tmp.evaluate(reduced);
			// move back to the original field
			res >>= *parent_field;
			return res;
		}
		// otherwise throw exception
		else throw IllegalCoercionException();
	}

}