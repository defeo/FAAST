#include "Types.hpp"
#include "Field.hpp"
#include <cstdlib>

using namespace std;
using namespace AS;

typedef Field<ZZ_p_Algebra> gfp;
typedef Field<zz_p_Algebra> GFp;
typedef Field<GF2_Algebra>  GFp2;
typedef FieldElement<ZZ_p_Algebra> gfp_E;
typedef FieldElement<zz_p_Algebra> GFp_E;
typedef FieldElement<GF2_Algebra>  GFp2_E;

int main(int argv, char* argc[]) {
	double cputime;

	cout << "Using " << gfp::Infrastructure::name << endl << endl;

	cputime = -NTL::GetTime();
	const gfp* K = &(gfp::createField(3,2));
	cputime += NTL::GetTime();
cout << gfp::Infrastructure::GFpE::modulus() << endl;
	cout << *K << " in " << cputime << endl;
	cout << "Time spent building the irreducible polynomial : "
		<< gfp::TIME.BUILDIRRED << endl << endl;
	
	for (int i = 1 ; i <= 4 ; i++) {
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << *K << " in " << cputime << endl;
cout << gfp::Infrastructure::GFpE::modulus() << endl;
	}
	
	cout << endl << "Time spent building the cyclotomic polynomial : "
		<< gfp::TIME.CYCLOTOMIC << endl << endl;
	
	for (int i = 1 ; i <= 10 ; i++) {
		gfp_E a = K->subField().random();
		vector<gfp_E> down;
		cputime = -GetTime();
		pushDown(a, down);
		cputime += GetTime();
	
		cout << a << endl;
		vector<gfp_E>::iterator it;
		for (it = down.begin() ; it != down.end(); it++)
			cout << *it << " ";
		cout << endl << "Pushdown computed in " << cputime << endl << endl;
	}
}
