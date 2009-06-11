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

typedef Field<ZZ_p_Algebra> GFp;
typedef Field<zz_p_Algebra> gfp;
typedef Field<GF2_Algebra>  GFp2;
typedef FieldElement<ZZ_p_Algebra> GFp_E;
typedef FieldElement<zz_p_Algebra> gfp_E;
typedef FieldElement<GF2_Algebra>  GFp2_E;

int main(int argv, char* argc[]) {
	double cputime;

	gfp::Infrastructure::BigInt p;
	long d, l, t;
	cin >> p; cin >> d; cin >> l; cin >> t;

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

	for (int i = 1 ; i <= l ; i++) {
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << *K << " in " << cputime << endl;

		for (int j = 0 ; j < i+t ; j++) {
			gfp_E a = K->random(), b;
			long n = (j>=i) ?
				d + RandomBnd(K->d - d) : d*power_long(p, j);
			double frobtime, pseudotime, naivetime;

			frobtime = -NTL::GetTime();
			b = a.frobenius(n);
			frobtime += NTL::GetTime();
			cout << "Fast " << n << "-th Frobenius in "
				<< frobtime << endl;

			pseudotime = -NTL::GetTime();
			b = a.pseudotrace(n);
			pseudotime += NTL::GetTime();
			cout << "Fast " << n << "-th Pseudotrace in "
				<< pseudotime << endl;

#ifdef FAAST_TIMINGS
			cout << "Pseudotraces precomputed in " <<
				gfp::TIME.PSEUDOTRACES << endl;
#endif

			naivetime = -NTL::GetTime();
			for (long i = 0 ; i < 10 ; i++)
				a.self_frobenius();
			naivetime += NTL::GetTime();
			double average = naivetime / 10;
			cout << "Naive Frobenius in " << average
				<< ", projected execution time " <<
				average * n << endl;

			cout << "Projected thresholds :" << endl
				<< "\tfrobenius " << floor(frobtime / average) << endl
				<< "\tpseudotrace " << floor(pseudotime / average) << endl;
		}
		cout << endl;
	}
}
