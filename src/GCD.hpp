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
		if (!P.parent_field && !Q.parent_field) return U = V = P;
		else if (!P.parent_field) {
			U = V = Q.parent_field->one();
			return Q;
		} else if (!Q.parent_field) {
			U = V = P.parent_field->one();
			return P;
		}
		P.sameLevel(Q);
		P.parent_field->switchContext();

		G.parent_field = U.parent_field = V.parent_field = P.parent_field;
		G.base = U.base = V.base = P.base;
		if (P.base) NTL::XGCD(G.repBase, U.repBase, V.repBase, P.repBase, Q.repBase);
		else NTL::XGCD(G.repExt, U.repExt, V.repExt, P.repExt, Q.repExt);
	}

	/* Half GCD between P and Q */
	template <class T> void
	HalfGCD(FieldPolynomial<T>& U0, FieldPolynomial<T>& V0,
			FieldPolynomial<T>& U1, FieldPolynomial<T>& V1,
			const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q)
	throw(NotInSameFieldException){
//TODO
	}
}
