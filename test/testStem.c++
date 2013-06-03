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
	\example testStem.c++
	This example illustrates how to use FAAST::Field::ArtinSchreierExtension(),
	pushDown() and liftUp().
*/

#include <faast.hpp>
#include <cstdlib>

using namespace std;
using namespace FAAST;

typedef Field<GF2_Algebra> gfp;
typedef FieldElement<GF2_Algebra> gfp_E;

int main(int argv, char* argc[]) {
	double cputime, totaltime;
	int retval = 0;

	gfp::Infrastructure::BigInt p;
	long d, l;
	if (cin.peek() != EOF) {
	  cin >> p; cin >> d; cin >> l;
	} else {
	  p = 2; d = 1; l = 4;
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

	cout << "\tCreate\tPushDow\tLiftUp\tPreLift\tMul\tInv" << endl;
	totaltime = -GetTime();
	for (int i = 1 ; i <= l ; i++) {
		cout << i << "\t";
		/** Construction **/
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		gfp_E a = K->random(), b;
		vector<gfp_E> down;

		/** Push-down / Lift-up **/
		cputime = -GetTime();
		pushDown(a, down);
		cputime += GetTime();
		cout << cputime << "\t";

		cputime = -GetTime();
		liftUp(down, b);
		cputime += GetTime();
		cout << cputime << "\t";
#ifdef FAAST_TIMINGS
		cout << gfp::TIME.LIFTUP;
#endif
		cout << "\t";

		if (a != b) {
			cout << "ERROR : Results don't match" << endl;
			cout << a << endl << b << endl;
			vector<gfp_E>::iterator it;
			for (it = down.begin() ; it != down.end() ; it++)
				cout << *it << " ";
			cout << endl;
			retval = 1;
		}

		/** Multiplication **/
		cputime = -GetTime();
		a*b;
		cputime += GetTime();
		cout << cputime << "\t";

		/** Inversion **/
		cputime = -GetTime();
		if (a != 0) a.inv();
		cputime += GetTime();
		cout << cputime << endl;
	}
	totaltime += GetTime();

#ifdef FAAST_TIMINGS
	cout << endl << "Time spent building the cyclotomic polynomial : "
		<< gfp::TIME.CYCLOTOMIC << endl;
#endif
	cout << endl;
	cout << "Total duration : " << totaltime << endl;

	return retval;
}
