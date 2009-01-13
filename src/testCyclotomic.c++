#include "Types.hpp"
#include "Field.hpp"
#include "utilities.hpp"
#include <cstdlib>

using namespace std;
using namespace AS;

typedef Field<ZZ_p_Algebra> GFp;
typedef Field<zz_p_Algebra> gfp;
typedef Field<GF2_Algebra>  GFp2;

int main(int argv, char* argc[]) {
	const GFp& K = GFp::createField(to_ZZ(3),1);
	const gfp& k = gfp::createField(5,1);
	const GFp2& L = GFp2::createField(2,1);

	while (true) {
		long n; cin >> n;
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
		cout << cputime << endl;

		GF2_Algebra::GFpX Phi2;
		cputime = -NTL::GetTime();
		cyclotomic<GF2_Algebra>(Phi2, n, L.characteristic());
		cputime += NTL::GetTime();
		//cout << Phi2 << endl;
		cout << cputime << endl;

		ZZ_p_Algebra::GFpX Phi3;
		cputime = -NTL::GetTime();
		cyclotomic<ZZ_p_Algebra>(Phi3, n, K.characteristic());
		cputime += NTL::GetTime();
		//cout << Phi3 << endl;
		cout << cputime << endl;
	}
}
