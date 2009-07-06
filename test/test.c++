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
	\example test.c++
	This example illustrates how to use FAAST::FieldElement::minimalPolynomials(),
	FAAST::FieldElement::affineMinimalPolynomial() and FAAST::FieldElement::evaluate().
*/

#include <faast.hpp>
#include <cstdlib>

using namespace std;
using namespace FAAST;

typedef Field<zz_p_Algebra>  gfp;
typedef FieldElement<zz_p_Algebra>  gfp_E;
typedef FieldPolynomial<zz_p_Algebra>  gfp_X;

int main(int argv, char* argc[]) {
	double cputime;
	int retval = 0;

	gfp::Infrastructure::BigInt p;
	long d, l;
	if (cin.peek() != EOF) {
	  cin >> p; cin >> d; cin >> l;
	} else {
	  p = 3; d = 2; l = 3;
	}

	cout << "Using " << gfp::Infrastructure::name << endl << endl;
	cputime = -NTL::GetTime();
	const gfp* K = &(gfp::createField(p,d));
	cputime += NTL::GetTime();
	cout << *K << " in " << cputime << endl;
#ifdef FAAST_TIMINGS
	cout << "Time spent building the irreducible polynomial : "
		<< gfp::TIME.BUILDIRRED << endl;
#endif
	cout << endl;

	cout << "\tCreate\tMinPol\tInterp\tEval" << endl;
	for (int i = 1 ; i <= l ; i++) {
		cout << i << "\t";
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		gfp_E a, b, c;
		a = K->random();
		b = K->random();

		vector<gfp_X> minpols;
		cputime = -NTL::GetTime();
		a.minimalPolynomials(K->baseField(),minpols);
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		gfp_X poly;
		bool nopol = false;
		cputime = -NTL::GetTime();
		try {
			poly = a.affineMinimalPolynomial(K->baseField(), b, minpols);
		} catch (NoSuchPolynomialException e) {
			nopol = true;
			cout << "*";
		}
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		if (!nopol) {
			cputime = -NTL::GetTime();
			c = a.evaluate(poly, minpols);
			cputime += NTL::GetTime();
			cout << cputime << endl;

			if (c != b) {
				cout << "ERROR 1 : Results don't match" << endl;
				cout << a << endl << b << endl << c << endl;
				cout << poly << endl;
				retval = 1;
			}

			poly >>= *K;
			c = poly.evaluate(a);
			if (c != b) {
				cout << "ERROR 2 : Results don't match" << endl;
				cout << a << endl << b << endl << c << endl;
				cout << poly << endl;
				retval = 1;
			}
		} else cout << "*\t" << endl;
	}
	return retval;
}
