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
/****************** GCD auxiliary routines ******************/
/* Half GCD between P and Q, assumes deg P >= deg Q */
template <class T> void
IterHalfGCD(FieldPolynomial<T>& U0, FieldPolynomial<T>& V0,
			FieldPolynomial<T>& U1, FieldPolynomial<T>& V1,
			FieldPolynomial<T>& P, FieldPolynomial<T>& Q,
			const long d) {

	U0 = V1 = P.parent().one();
	U1 = V0 = P.parent().zero();

	long goal = P.degree() - d;

	if (Q.degree() <= goal)
		return;

	FieldPolynomial<T> q, t;

	while (Q.degree() > goal) {
		q.division(P, Q); t.mod(P, Q);
		P = Q; Q = t;

		t.product(q, U1); t.negate(); t += U0;
		U0 = U1; U1 = t;

		t.product(q, V1); t.negate(); t += V0;
		V0 = V1; V1 = t;
	}
}



/* Half GCD between P and Q, assumes deg P > deg Q */
template <class T> void
RecHalfGCD(FieldPolynomial<T>& U0, FieldPolynomial<T>& V0,
		FieldPolynomial<T>& U1, FieldPolynomial<T>& V1,
		const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q,
		const long d) {

	if (Q == 0 || Q.degree() <= P.degree() - d) {
		U0 = V1 = P.parent().one();
		U1 = V0 = P.parent().zero();

		return;
	}

	long n = P.degree() - 2*d + 2;
	if (n < 0) n = 0;

	FieldPolynomial<T> P1, Q1;
	P1.RightShift(P, n);
	Q1.RightShift(Q, n);

	if (d*P.parent().degree() <= T::consts.HalfGCD_CROSSOVER) {
		IterHalfGCD<T>(U0, V0, U1, V1, P1, Q1, d);
		return;
	}
	long d1 = (d + 1)/2;
	if (d1 >= d) d1 = d - 1;

	FieldPolynomial<T> u0, v0, u1, v1;

	RecHalfGCD<T>(u0, v0, u1, v1, P1, Q1, d1);
	// matrix-vector product
	// |P1|    |u0 v0| |P1|
	// |Q1| <- |u1 v1| |Q1|
	FieldPolynomial<T> t1, t2;
	t1.product(u0, P1);
	t2.product(v0, Q1);
	t1.sum(t1, t2);
	t2.product(u1, P1);
	Q1.product(v1, Q1);
	Q1.sum(t2, Q1);
	P1 = t1;

	long d2 = Q1.degree() - P.degree() + n + d;

	if (Q1 == 0 || d2 <= 0) {
		U0 = u0; V0 = v0; U1 = u1; V1 = v1;
		return;
	}

	FieldPolynomial<T> q, y0, w0, y1, w1;

	q.division(P1, Q1); P1 %= Q1;
	RecHalfGCD<T>(y0, w0, y1, w1, Q1, P1, d2);

	// matrix multiplication
	// |u0 v0|    |0  1| |u0 v0|
	// |u1 v1| <- |1 -q| |u1 v1|
	t1.product(q, u1); t1.negate(); t1 += u0;
	u0 = u1; u1 = t1;
	t1.product(q, v1); t1.negate(); t1 += v0;
	v0 = v1; v1 = t1;

	// matrix multiplication (Strassen formula)
	// |U0 V0|   |y0 w0| |u0 v0|
	// |U1 V1| = |y1 w1| |u1 v1|
	U0 = V0 = U1 = V1 = P.parent().zero();

	t1.difference(y0, w0); t1 *= v1;
	U0 += t1; V0 -= t1;

	t1.difference(y1, w1); t1 *= u0;
	U1 += t1; V1 -= t1;

	t1.sum(u0, u1); t1 *= w1;
	U0 -= t1; U1 += t1;

	t1.sum(v0, v1); t1 *= y0;
	V0 += t1; V1 -= t1;

	t2.sum(y0, w0); t1.difference(v1, u0); t1 *= t2;
	U0 -= t1; V1 += t1;

	t2.sum(y0, y1); t1.sum(u0, v0); t1 *= t2;
	V1 += t1;

	t2.sum(w0, w1); t1.sum(u1, v1); t1 *= t2;
	U0 += t1;
}



/****************** GCD ******************/
	/* GCD between P and Q */
	template <class T> FieldPolynomial<T>
	GCD(const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q)
	throw(NotInSameFieldException) {
		if (!P.parent_field) return Q;
		if (!Q.parent_field) return P;
		P.sameLevel(Q);
		P.parent_field->switchContext();

		FieldPolynomial<T> res = P.parent_field->zero();
		if (P.base) NTL::GCD(res.repBase, P.repBase, Q.repBase);
		else {
			FieldPolynomial<T> U, V;
			XGCD<T>(P, Q, U, V, res);
		}

		return res;
	}

	/* XGCD between P and Q */
	template <class T> void
	XGCD(const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q,
	FieldPolynomial<T>& U, FieldPolynomial<T>& V, FieldPolynomial<T>& G)
	throw(NotInSameFieldException) {
		if (!P.parent_field && !Q.parent_field) G = U = V = P;
		else if (!P.parent_field) {
			U = Q.parent_field->zero();
			V = Q.parent_field->one();
			G = Q;
		} else if (!Q.parent_field) {
			U = P.parent_field->one();
			V = P.parent_field->zero();
			G = P;
		} else {
			P.sameLevel(Q);
			FieldPolynomial<T> U1, V1;
			if (P.degree() >= Q.degree()) {
				RecHalfGCD<T>(U, V, U1, V1, P, Q, P.degree()+1);
			} else {
				RecHalfGCD<T>(V, U, V1, U1, Q, P, Q.degree()+1);
			}

			G = U*P + V*Q;
			// normalize
			FieldElement<T> lc; G.getCoeff(G.degree(), lc);
			G /= lc; U /= lc; V /= lc;
		}
	}

	/* Half GCD between P and Q */
	template <class T> void
	HalfGCD(FieldPolynomial<T>& U0, FieldPolynomial<T>& V0,
			FieldPolynomial<T>& U1, FieldPolynomial<T>& V1,
			const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q,
			const long d)
	throw(NotInSameFieldException, BadParametersException){
		if (d > P.degree() && d > Q.degree())
			throw BadParametersException("Negative goal degree for HalfGCD");
		if (d < 0)
			throw BadParametersException("Negative parameter for HalfGCD");

		if (!P.parent_field && !Q.parent_field) U0 = V0 = U1 = V1 = P;
		else if (!P.parent_field) {
			U0 = V1 = Q.parent_field->zero();
			U1 = V0 = Q.parent_field->one();
		} else if (!Q.parent_field) {
			U1 = V0 = Q.parent_field->zero();
			U0 = V1 = P.parent_field->one();
		} else {
			P.sameLevel(Q);

			if (P.degree() >= Q.degree()) {
				RecHalfGCD<T>(U0, V0, U1, V1, P, Q, d);
			} else {
				RecHalfGCD<T>(V0, U0, V1, U1, Q, P, d);
			}
		}
	}
}
