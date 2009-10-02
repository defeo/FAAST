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

NTL_OPEN_NNS
//static void HalfGCD(_NTL_GF2XMatrix& M_out, const GF2X& U, const GF2X& V, long d_red);
NTL_CLOSE_NNS


namespace FAAST {
/****************** GCD auxiliary routines ******************/
/* Half GCD between P and Q, assumes deg P >= deg Q */
template <class T> void
IterHalfGCD(typename T::GFpEX& U0, typename T::GFpEX& V0,
			typename T::GFpEX& U1, typename T::GFpEX& V1,
			const typename T::GFpEX& P, const typename T::GFpEX& Q,
			const long d) {
	typedef typename T::GFpEX   GFpEX;


	U0.SetMaxLength(d);
	V0.SetMaxLength(d);
	U1.SetMaxLength(d);
	V1.SetMaxLength(d);

	set(U0); clear(V0);
	clear(U1); set(V1);

	long goal = deg(P) - d;

	if (deg(Q) <= goal)
		return;

	GFpEX q, t;
	t.SetMaxLength(d);

	while (deg(Q) > goal) {
		PlainDivRem(q, P, P, Q);
		swap(P, Q);

		mul(t, q, U1);
		sub(t, U0, t);
		U0 = U1;
		U1 = t;

		mul(t, q, V1);
		sub(t, V0, t);
		V0 = V1;
		V1 = t;
	}
}



/* Half GCD between P and Q, assumes deg P > deg Q */
template <class T> void
HalfGCD(typename T::GFpEX& U0, typename T::GFpEX& V0,
		typename T::GFpEX& U1, typename T::GFpEX& V1,
		const typename T::GFpEX& P, const typename T::GFpEX& Q,
		const long d) {
	typedef typename T::GFpE    GFpE;
	typedef typename T::GFpEX   GFpEX;

	if (IsZero(Q) || deg(Q) <= deg(P) - d) {
		set(U0); clear(V0);
		clear(U1); set(V1);

		return;
	}

	long n = deg(P) - 2*d + 2;
	if (n < 0) n = 0;

	GFpEX P1, Q1;

	RightShift(P1, P, n);
	RightShift(Q1, Q, n);

	if (d*GFpE::degree()*2 <= NTL_zz_pX_HalfGCD_CROSSOVER) {
		IterHalfGCD(U0, V0, U1, V1, P1, Q1, d);
		return;
	}
	long d1 = (d+ 1)/2;
	if (d1 < 1) d1 = 1;
	if (d1 >= d) d1 = d - 1;

	GFpEX u0, v0, u1, v1;

	HalfGCD(u0, v0, u1, v1, P1, Q1, d1);
	// matrix multiplication
	mul(P1, Q1, M1);

	long d2 = deg(Q1) - deg(P) + n + d;

	if (IsZero(Q1) || d2 <= 0) {
		U0 = u0; V0 = v0; U1 = u1; V1 = v1;
		return;
	}

	GFpEX q;
	GFpEX y0, w0, y1, w1;

	DivRem(q, P1, P1, Q1);
	swap(P1, Q1);
	HalfGCD(y0, w0, y1, w1, P1, Q1, d2);

	GFpEX t(INIT_SIZE, deg(u1)+deg(q)+1);

	mul(t, q, u1);
	sub(t, u0, t);
	swap(u0, u1);
	swap(u1, t);

	t.kill();

	t.SetMaxLength(deg(v1)+deg(q)+1);

	mul(t, q, v1);
	sub(t, v0, t);
	swap(v0, v1);
	swap(v1, t);

	t.kill();

	// matrix multiplication
	mul(M_out, M2, M1);
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
		else NTL::GCD(res.repExt, P.repExt, Q.repExt);

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
			V = Q.parent_field->zero();
			G = P;
		} else {
			P.sameLevel(Q);
			P.parent_field->switchContext();

			if (P.base) {
				NTL::XGCD(G.repBase, U.repBase, V.repBase, P.repBase, Q.repBase);
				G.repExt = U.repExt = V.repExt = 0;
			} else {
				NTL::XGCD(G.repExt, U.repExt, V.repExt, P.repExt, Q.repExt);
				G.repBase = U.repBase = V.repBase = 0;
			}
			G.parent_field = U.parent_field = V.parent_field = P.parent_field;
			G.base = U.base = V.base = P.base;
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

		if (!P.parent_field && !Q.parent_field) return U0 = V0 = U1 = V1 = P;
		else if (!P.parent_field) {
			U0 = V1 = Q.parent_field->zero();
			U1 = V0 = Q.parent_field->one();
		} else if (!Q.parent_field) {
			U1 = V0 = Q.parent_field->zero();
			U0 = V1 = P.parent_field->one();
		} else {
			P.sameLevel(Q);
			P.parent_field->switchContext();

			if (P.base) {
				typename T::GFpX4Mat M;
				if (P.degree() >= Q.degree()) {
					NTL::HalfGCD(M, P.repBase, Q.repBase, d);
					U0.repBase = M(0,0); V0.repBase = M(0,1);
					U1.repBase = M(1,0); V1.repBase = M(1,1);
				} else {
					NTL::HalfGCD(M, Q.repBase, P.repBase, d);
					U0.repBase = M(0,1); V0.repBase = M(0,0);
					U1.repBase = M(1,1); V1.repBase = M(1,0);
				}
				U0.repExt = V0.repExt = U1.repExt = V1.repExt = 0;
			} else {
				if (P.degree() >= Q.degree()) {
					HalfGCD(U0.repExt, V0.repExt, U1.repExt, V1.repExt, P.repExt, Q.repExt, d);
				} else {
					HalfGCD(V0.repExt, U0.repExt, V1.repExt, U1.repExt, Q.repExt, P.repExt, d);
				}
				U0.repBase = V0.repBase = U1.repBase = V1.repBase = 0;
			}
			U0.base = V0.base = U1.base = V1.base = P.base;
			U0.parent_field = V0.parent_field = U1.parent_field = V1.parent_field = P.parent_field;
		}
	}
}
