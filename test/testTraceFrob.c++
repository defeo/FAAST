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
	\example testTraceFrob.c++
	This example illustrates how to use FAAST::FieldElement::frobenius() and
	FAAST::FieldElement::pseudotrace().
*/

#include <faast.hpp>
#include <cstdlib>

using namespace std;
using namespace FAAST;

typedef Field<zz_p_Algebra> gfp;
typedef FieldElement<zz_p_Algebra> gfp_E;

int main(int argv, char* argc[]) {
	double cputime;

	gfp::Infrastructure::BigInt p;
	long d, l, t;
	if (cin.peek() != EOF) {
	  cin >> p; cin >> d; cin >> l; cin >> t;
	} else {
	  p = 2; d = 3; l = 4; t = 0;
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

	cout << "\t\tFrob\tPTr\tPrePTr\tNFrob\tNProj" << endl;
	for (int i = 1 ; i <= l ; i++) {
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << i << "\t" << cputime << endl;

		for (int j = 0 ; j < i+t ; j++) {
			gfp_E a = K->random(), b;
			long n = (j>=i) ?
				d + RandomBnd(K->d - d) : d*power_long(p, j);
			double frobtime, pseudotime, naivetime;

			frobtime = -NTL::GetTime();
			b = a.frobenius(n);
			frobtime += NTL::GetTime();
			cout << "\t" << n << "\t" 
				<< frobtime << "\t";

			pseudotime = -NTL::GetTime();
			b = a.pseudotrace(n);
			pseudotime += NTL::GetTime();
			cout << pseudotime << "\t";

#ifdef FAAST_TIMINGS
			cout << gfp::TIME.PSEUDOTRACES << "\t";
#endif

			naivetime = -NTL::GetTime();
			for (long i = 0 ; i < 10 ; i++)
				a.self_frobenius();
			naivetime += NTL::GetTime();
			double average = naivetime / 10;
			cout <<  average << "\t" <<
				average * n << endl;
		}
	}
	return 0;
}
