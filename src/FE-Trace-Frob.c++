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
#import "FAAST/utilities.hpp"
#ifdef FAAST_DEBUG
#import <string>
#import <sstream>
#endif

namespace FAAST {
/****************** Arithmetics ******************/
	/* n-th iterated frobenius. This algorithm is detailed in the
	 * long version.
	 */
	template <class T> void
	FieldElement<T>::self_frobenius(long n) throw() {
		if (isZero()) return;
		if (base) return;

		n %= parent_field->d;
		if (n < 0) n = parent_field->d + n;
		if (n == 0) return;
		long p = parent_field->p;
		long smalld = parent_field->baseField().d;
		// The small part
		const long r = n % smalld;
		SmallFrob(r);
		// The big part
		n /= smalld;
		long j = 0;
		while (n != 0) {
			const long c = n % p;
			for (long i = 0 ; i < c ; i++) BigFrob(j);
			n /= p ; j++;
		}
	}

	template <class T> void
	FieldElement<T>::self_trace(const Field<T>& F)
	throw(NotASubFieldException) {
		// zero stays zero
		if (!parent_field) {
			*this = F.zero();
		}
		// handle errors
		else if (!F.isSubFieldOf(*parent_field)) {
			throw NotASubFieldException();
		}
		// if F is a prime field, this is the absolute trace
		else if (F.isPrimeField()) {
			self_trace();
		}
		// else, the trace is minus the coefficient of x_i^(p-1)
		else {
			// move to the stem
			parent_field = parent_field->stem;

			vector<FieldElement<T> > down;
			while (parent_field != F.stem) {
				FAAST::pushDown(*this, down);
				if (down.size() == unsigned(parent_field->p)) {
					*this = down[long(parent_field->p) - 1];
					negate();
				} else *this = F.stem->zero();
			}
			// move out of the stem
			parent_field = &F;
		}
	}

	template <class T> void
	FieldElement<T>::self_pseudotrace(unsigned long n) throw() {
		if (isZero() || n == 1) return;
		if (n == 0) *this = parent_field->zero();
		if (base) {
			*this *= parent_field->scalar(n);
			return;
		}

		// first remove the part above d
		FieldElement<T> traces;
		if (n >= unsigned(parent_field->d)) {
			traces = trace();
			traces *= parent_field->scalar(n / parent_field->d);
		}
		// now go under d
		n %= parent_field->d;

		long p = parent_field->p;
		long smalld = parent_field->baseField().d;
		// The small part
		const long r = n % smalld;
		// The big part
		n /= smalld;
		const long j = NumPits(p, n);
		if (j > 0) {
			vector<FieldElement<T> > v;
			BigPTraceVector(v, j-1);
			long i;
			for (i = j-1 ; i >= 0 && n != 0 ; i--) {
				const long c = n / power_long(p, i);
				if (c == 0) {
					v[i] = 0;
				} else if (c > 1) {
					FieldElement<T> t = v[i];
					for (long h = 1 ; h < c ; h++) {
						t.BigFrob(i);
						v[i] += t;
					}
				}
				if (i < j-1) {
					for (long h = 0 ; h < c ; h++)
						v[i+1].BigFrob(i);
					v[i] += v[i+1];
				}
				n %= power_long(p, i);
			}
			v[i+1].SmallFrob(r);
			this->SmallPTrace(r);
			*this += v[i+1];
		} else {
			SmallPTrace(r);
		}
		// add the part above d, if needed
		*this += traces;
	}

/****************** Helpers for frobenius and trace ******************/
	/* See Section 5 and the long version (probably section 6). */

	/* p^j-th iterated frobenius */
	template <class T> void FieldElement<T>::BigFrob(const long j) {
#ifdef FAAST_DEBUG
		if (j < 0 || j >= parent_field->height) {
			stringstream msg;
			msg << "Bad input to BigFrob : " << j << ".";
			throw FAASTException(msg.str().c_str());
		}
#endif
		if (isScalar()) return;

		BigInt p = parent_field->p;
		// step 2
		vector<FieldElement<T> > down;
		FAAST::pushDown(*this, down);
		down.resize(p);
		// step 3
		if (j < parent_field->height - 1) {
			for (BigInt i = 0 ; i < p ; i++)
			down[i].BigFrob(j);
		}
		// step 5
		const FieldElement<T>& beta = parent_field->getPseudotrace(j);
		vector<FieldElement<T> > result;
		result.resize(p);
		for (BigInt i = 0 ; i < p ; i++) {
			result[i] = down[long(p)-1];
			for (long j = long(i) - 1 ; j >= 0 ; j--) {
				result[j] *= beta;
				if (j > 0) result[j] += result[j-1];
				else result[j] += down[long(p)-long(i)-1];
			}
		}
		// step 6
		const Field<T>* parent = parent_field;
		FAAST::liftUp(result, *this);
		parent_field = parent;
	}

	/* n-th iterated frobenius, n < d */
	template <class T> void FieldElement<T>::SmallFrob(const long n) {
#ifdef FAAST_DEBUG
		if (n < 0 || n >= parent_field->baseField().d) {
			stringstream msg;
			msg << "Bad input to SmallFrob : " << n << ".";
			throw FAASTException(msg.str().c_str());
		}
#endif
		if (isScalar()) return;

		for (long i = 0 ; i < n ; i++)
			self_frobenius();
	}

	/* p^j-th pseudotrace */
	template <class T> void FieldElement<T>::BigPTrace(const long j) {
#ifdef FAAST_DEBUG
		if (j < 0 || j >= parent_field->height) {
			stringstream msg;
			msg << "Bad input to BigPTrace : " << j << ".";
			throw FAASTException(msg.str().c_str());
		}
#endif
		if (isZero()) return;

		SmallPTrace(parent_field->baseField().d);
		for (long i = 1 ; i <= j ; i++) {
			FieldElement<T> t = *this;
			for (BigInt h = 1 ; h < parent_field->p ; h++) {
				t.BigFrob(i-1);
				*this += t;
			}
		}
	}

	/* Put in the vector v all the p^id pseudotraces for 0 <= i <= j */
	template <class T> void
	FieldElement<T>::BigPTraceVector(vector<FieldElement<T> >& v,
	const long j) const {
#ifdef FAAST_DEBUG
		if (j < 0 || j >= parent_field->height) {
			stringstream msg;
			msg << "Bad input to BigPTraceVector : " << j << ".";
			throw FAASTException(msg.str().c_str());
		}
#endif
		v.clear(); v.resize(j+1);
		v[0] = *this;
		v[0].SmallPTrace(parent_field->baseField().d);
		for (long i = 1 ; i <= j ; i++) {
			v[i] = v[i-1];
			FieldElement<T> t = v[i];
			for (BigInt h = 1 ; h < parent_field->p ; h++) {
				t.BigFrob(i-1);
				v[i] += t;
			}
		}
	}

	/* n-th pseudotrace, n < d */
	template <class T> void FieldElement<T>::SmallPTrace(const long n) {
#ifdef FAAST_DEBUG
		if (n < 0 || n > parent_field->baseField().d) {
			stringstream msg;
			msg << "Bad input to SmallPTrace : " << n << ".";
			throw FAASTException(msg.str().c_str());
		}
#endif
		if (isZero()) return;

		if (n == 0) {
			*this = parent_field->zero();
		} else {
			FieldElement<T> t(*this);
			for (long i = 1 ; i < n ; i++) {
				t.self_frobenius();
				*this += t;
			}
		}
	}

}
