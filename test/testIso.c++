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
	\example testIso.c++
	This example illustrates how to navigate a \ref Field_lattices "field lattice"
	and how to use FAAST::Field::ArtinSchreierExtension(),
	FAAST::Field::toBivariate() and FAAST::Field::toUnivariate().

*/

#include <faast.hpp>
#include <cstdlib>

using namespace std;
using namespace FAAST;

typedef Field<ZZ_p_Algebra>  gfp;
typedef FieldElement<ZZ_p_Algebra>  gfp_E;

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

	cout << "\tCreate\tCrStem\tPushDow\tLiftUp\tPrePseu\tPreLift" << endl;
	for (int i = 1 ; i <= l ; i++) {
		cout << i << "\t";
		gfp_E alpha;
		do {
			alpha = K->random();
		} while (alpha.trace() == 0);
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension(alpha));
		cputime += NTL::GetTime();
		cout << cputime << "\t";
#ifdef FAAST_TIMINGS
		cout << gfp::TIME.BUILDSTEM;
#endif
		cout << "\t";

		const gfp& L = K->stemField().subField();
		vector<gfp_E> v;
		gfp_E a = K->random(), b;
		cputime = -NTL::GetTime();
		L.toBivariate(a, v);
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		cputime = -NTL::GetTime();
		K->toUnivariate(v, b);
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		if (a != b) {
			cout << "ERROR : Results don't match" << endl;
			cout << a << endl;
			cout << b << endl;
			vector<gfp_E>::iterator it;
			for (it = v.begin() ; it != v.end() ; it++)
				cout << *it << " ";
			cout << endl;
			retval = 1;
		}

#ifdef FAAST_TIMINGS
		cout << gfp::TIME.PSEUDOTRACES << "\t";
		cout << gfp::TIME.LIFTUP;
#endif
		cout << endl;
	}
#ifdef FAAST_TIMINGS
	cout << endl << "Time spent inverting the matrix " <<
		gfp::TIME.ARTINMATRIX << endl;
#endif

	return retval;
}
