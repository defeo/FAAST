#include "Types.hpp"
#include "Field.hpp"
#include <cstdlib>

using namespace std;
using namespace AS;

typedef Field<ZZ_p_Algebra> GFp;
typedef Field<zz_p_Algebra> GFp2;
typedef Field<GF2_Algebra>  gfp;
typedef FieldElement<ZZ_p_Algebra> GFp_E;
typedef FieldElement<zz_p_Algebra> GFp2_E;
typedef FieldElement<GF2_Algebra>  gfp_E;

int main(int argv, char* argc[]) {
	double cputime;
	
	gfp::Infrastructure::BigInt p;
	long d, l;
	cin >> p; cin >> d; cin >> l;

	cout << "Using " << gfp::Infrastructure::name << endl << endl;
	cputime = -NTL::GetTime();
	const gfp* K = &(gfp::createField(p,d));
	cputime += NTL::GetTime();
	cout << *K << " in " << cputime << endl;
	cout << "Time spent building the irreducible polynomial : "
		<< gfp::TIME.BUILDIRRED << endl << endl;
	
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
		cout << gfp::TIME.BUILDSTEM << "\t";
			
		const gfp& L = K->stemField().subField();
		vector<gfp_E> v;
		gfp_E a = K->random(), b;
		cputime = -NTL::GetTime();
		L.pushDown(a, v);
		cputime += NTL::GetTime();
		cout << cputime << "\t";
		
		cputime = -NTL::GetTime();
		K->liftUp(v, b);
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
		}
		
		cout << gfp::TIME.PSEUDOTRACES << "\t";
		cout << gfp::TIME.LIFTUP << endl;
	}
	cout << endl << "Time spent inverting the matrix" <<
		gfp::TIME.ARTINMATRIX << endl;
}
