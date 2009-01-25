#include "Types.hpp"
#include "Field.hpp"
#include <cstdlib>

using namespace std;
using namespace AS;

typedef Field<ZZ_p_Algebra> GFp;
typedef Field<zz_p_Algebra> gfp;
typedef Field<GF2_Algebra>  GFp2;
typedef FieldElement<ZZ_p_Algebra> GFp_E;
typedef FieldElement<zz_p_Algebra> gfp_E;
typedef FieldElement<GF2_Algebra>  GFp2_E;

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
	
	for (int i = 1 ; i <= l ; i++) {
		gfp_E alpha;
		do {
			alpha = K->random();
		} while (alpha.trace() == 0);
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension(alpha));
		cputime += NTL::GetTime();
		cout << *K << " in " << cputime << endl;
		cout << "Time spent building the stem " <<
			gfp::TIME.BUILDSTEM << endl;
			
		const gfp& L = K->stemField().subField();
		vector<gfp_E> v;
		gfp_E a = K->random(), b;
		cputime = -NTL::GetTime();
		L.pushDown(a, v);
		cputime += NTL::GetTime();
		cout << "Push-down in " << cputime << endl;
		
		cputime = -NTL::GetTime();
		K->liftUp(v, b);
		cputime += NTL::GetTime();
		cout << "Lift-up in " << cputime << endl;
		
		if (a != b) {
			cout << "ERROR : Results don't match" << endl;
			cout << a << endl;
			cout << b << endl;
			vector<gfp_E>::iterator it;
			for (it = v.begin() ; it != v.end() ; it++)
				cout << *it << " ";
			cout << endl;
		}
		
		cout << "Time spent in pseudotrace precomputation " <<
			gfp::TIME.PSEUDOTRACES << endl;
		cout << "Time spent in Lift-up precomputation " <<
			gfp::TIME.LIFTUP << endl << endl;
	}
	cout << "Time spent inverting the matrix" <<
		gfp::TIME.ARTINMATRIX << endl;
}
