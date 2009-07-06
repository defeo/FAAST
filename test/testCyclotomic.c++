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
#include <FAAST/utilities.hpp>
#include <cstdlib>

using namespace std;
using namespace FAAST;

typedef Field<ZZ_p_Algebra> GFp;
typedef Field<zz_p_Algebra> gfp;
typedef Field<GF2_Algebra>  GFp2;

int main(int argv, char* argc[]) {
	const GFp& K = GFp::createField(to_ZZ(3),1);
	const gfp& k = gfp::createField(5,1);
	const GFp2& L = GFp2::createField(2,1);

	long n;
	if (cin.peek() != EOF) cin >> n;
	else n = 10;

	vector<pair<long,int> > factors;
	factor(n, factors);
	vector<pair<long,int> >::iterator it;
	for (it = factors.begin() ; it != factors.end() ; it++) {
	  cout << "(" << it->first << ", " << it->second << ") ";
	}
	cout << endl;
	zz_p_Algebra::GFpX Phi;
	double cputime = -NTL::GetTime();
	cyclotomic<zz_p_Algebra>(Phi, n, k.characteristic());
	cputime += NTL::GetTime();
	//cout << Phi << endl;
	cout << cputime << "\t";
	
	GF2_Algebra::GFpX Phi2;
	cputime = -NTL::GetTime();
	cyclotomic<GF2_Algebra>(Phi2, n, L.characteristic());
	cputime += NTL::GetTime();
	//cout << Phi2 << endl;
	cout << cputime << "\t";
	
	ZZ_p_Algebra::GFpX Phi3;
	cputime = -NTL::GetTime();
	cyclotomic<ZZ_p_Algebra>(Phi3, n, K.characteristic());
	cputime += NTL::GetTime();
	//cout << Phi3 << endl;
	cout << cputime << endl;

	return 0;
}
