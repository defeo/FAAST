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
	const GFp& K = GFp::createField(to_ZZ(5),25);
	const gfp& k = gfp::createField(3,5);
	const GFp2& L = GFp2::createField(2,3);
	const GFp& base = K.baseField();

/*	cout << K.characteristic() << " " << K.degree() << endl 
		<< k.characteristic() << " " << k.degree() << endl
		<< L << ", " << base << endl
		<< L.cardinality() << " " << K.cardinality() << " "
		<< base.cardinality() << " " << k.cardinality() << endl;
	cout << K.isOverFieldOf(base) << K.isIsomorphic(K) << (K==K) << endl;*/

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
	}
}
