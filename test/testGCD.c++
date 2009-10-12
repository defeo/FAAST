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
#include <faast.hpp>
#include <cstdlib>

namespace FAAST {
	template <class T> void
	IterHalfGCD(FieldPolynomial<T>& U0, FieldPolynomial<T>& V0,
				FieldPolynomial<T>& U1, FieldPolynomial<T>& V1,
				FieldPolynomial<T>& P, FieldPolynomial<T>& Q,
				const long d);
	template <class T> void
	RecHalfGCD(FieldPolynomial<T>& U0, FieldPolynomial<T>& V0,
				FieldPolynomial<T>& U1, FieldPolynomial<T>& V1,
				const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q,
				const long d);
}

using namespace std;
using namespace FAAST;

typedef Field<GF2_Algebra> GFp;
typedef FieldPolynomial<GF2_Algebra> GFpX;

int main(int argv, char* argc[]) {
	double cputime;
	long p, n, d;
	if (cin.peek() != EOF) {
		cin >> p; cin >> n; cin >> d;
	} else {
		p = 2; n = 10; d = 10;
	}
	
	const GFp& K = GFp::createField(p, n);

	GFpX P, Q, P1, Q1, G1, G2, U0, V0, U1, V1;
	for (long i = 0 ; i < d ; i++) {
		P.setCoeff(i, K.random());
		Q.setCoeff(i, K.random());
	}
	P.setCoeff(d);

	for (long t = 0 ; t <= d ; t += 10) {
		GFp::Infrastructure::consts.HalfGCD_CROSSOVER = t;

		cputime = -NTL::GetTime();
		RecHalfGCD(U0, V0, U1, V1, P, Q, d/2);
		G1 = U0*P + V0*Q;
		cputime += NTL::GetTime();
		cout << t << "\t" << cputime << endl;
		cout.flush();
	}
/*
	P1 = P; Q1 = Q;
	cputime = -NTL::GetTime();
	IterHalfGCD(U0, V0, U1, V1, P1, Q1, d+1);
	G2 = U0*P + V0*Q;
	cputime += NTL::GetTime();
	cout << cputime << endl;
	cout.flush();
	
	if (G1 != G2) {
		cout << "ERROR : Results don't match" << endl;
		cout << P << endl << Q << endl << G1 << endl;
		cout << G2 << endl;
		return 1;
	} else return 0;*/
}
