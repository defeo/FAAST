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
	/* The algorithm ApproximateAS from Section 6 */
	template <class T> void approximateAS(
	vector<FieldElement<T> >& v, const FieldElement<T>& tr) {
		const Field<T>& parent = v[0].parent();
		for (long j = long(parent.p) - 1 ; j >= 1 ; j--) {
			FieldElement<T> t = tr;
			FieldElement<T> binom =
				parent.primeField().scalar(j);
			v[j] = v[j-1];
			for (long h = j + 1 ; h <= long(parent.p) - 1 ; h++) {
				t *= tr;
				binom *= parent.primeField().scalar(h);
				binom /= parent.primeField().scalar(h - j + 1);
				v[j] -= v[h] * (t * binom);
			}
			v[j] /= tr * parent.primeField().scalar(j);
			v[j] >>= parent;
		}
		v[0] = parent.zero();
	}

	/* Couveignes' algorithm (Section 6).
	 * Assumes Tr(alpha) = 0
	 */
	template <class T> void Field<T>::couveignes00(
	FieldElement<T>& res, const FieldElement<T>& alpha) const {
#ifdef FAAST_DEBUG
		if (alpha.trace() != 0)
			throw FAASTException("Bad input to couveignes00.");
#endif
		// step 0
		if (alpha.isZero()) {
			res = alpha;
			return;
		}

		const Field<T>& parent = *(alpha.parent_field);
		long i = parent.height;
		// step 1
		if (i == 0) {
			parent.switchContext();
			const MatGFp& artin = parent.getArtinMatrix();
#ifdef FAAST_DEBUG
			if (parent.stem->artinLine < 0)
				throw FAASTException("Bad Artin Matrix.");
#endif
			VecGFp low, high;
			VectorCopy(low, rep(alpha.repExt), parent.stem->artinLine);
			VectorCopy(high, RightShift(rep(alpha.repExt),
				parent.stem->artinLine+1),
				parent.d - 1 - parent.stem->artinLine);
			append(low, high);
			// apply the artin matrix
			VecGFp resV = artin * low;
			GFpX resX; conv(resX, resV);
			resX <<= 1;
			conv(res.repExt, resX);
			res.repBase = 0;
			res.base = false;
			res.parent_field = &parent;
			return;
		}
		// step 2
		FieldElement<T> eta = alpha;
		eta.BigPTrace(i-1);
		// step 3
		vector<FieldElement<T> > etas;
		FAAST::pushDown(eta, etas);
		etas.resize(p);
#ifdef FAAST_DEBUG
		if (etas[long(parent.p)-1] != 0)
			throw FAASTException("Error in couveignes00.");
#endif
		approximateAS<T>(etas, parent.stem->alpha->trace());
		// step 4
		vector<FieldElement<T> > alphas;
		FAAST::pushDown(alpha, alphas);
		alphas.resize(1);
		FieldElement<T> alpha1 = alphas[0];
		FieldElement<T> gamma = alphas[0].parent_field->one();
		for (BigInt j = 1 ; j < parent.p ; j++) {
			gamma *= *(parent.stem->alpha);
			FieldElement<T> tmp = etas[long(j)];
			tmp.self_frobenius();
			tmp *= gamma;
			alpha1 -= tmp;
		}
		// step5
		couveignes00(etas[0], alpha1);
		// step 6
		FAAST::liftUp(etas, res);
	}

	/* Return a root in this field of the polynomial
	 * 			X^p - X - alpha
	 * throws IllegalCoercionException if alpha can't be
	 * 		coerced to this field
	 * throws IsIrreducibleException if the polynomial has
	 * 		no roots in this field
	 */
	template <class T> FieldElement<T>
	Field<T>::Couveignes2000(const FieldElement<T>& alpha)
	const throw(IllegalCoercionException, IsIrreducibleException) {
		if (!isOverFieldOf(*(alpha.parent_field)))
			throw IllegalCoercionException();

		// if alpha lies in this level
		if (alpha.parent_field->stem == stem) {
			if (alpha.trace() != 0)
				throw IsIrreducibleException();
			FieldElement<T> res;
			couveignes00(res, alpha);
			return res >> *this;
		}
		// if alpha lies one level below
		else if (alpha.parent_field->stem->overfield == stem) {
			FieldElement<T> tr = alpha.trace();
			if (tr == 0) {
				FieldElement<T> res =
					stem->subfield->Couveignes2000(alpha);
				return res >> *this;
			} else {
				// simplified version of couveignes00
				vector<FieldElement<T> > v;
				v.resize(2);
				v[1] = tr / stem->alpha->trace();
				v[0] = *(stem->alpha);
				v[0] >>= *(alpha.parent_field);
				v[0] *= -v[1];
				v[0] += alpha;
				couveignes00(v[0], v[0]);
				v[1] >>= *(v[0].parent_field);
				FieldElement<T> res;
				FAAST::liftUp(v, res);
				return res >> *this;
			}
		}
		// else
		else {
			FieldElement<T> res =
				stem->subfield->Couveignes2000(alpha);
			return res >> *this;
		}
	}
}
