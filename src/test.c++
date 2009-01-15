#include "Types.hpp"
#include "Field.hpp"
#include <cstdlib>

using namespace std;
using namespace AS;

typedef Field<ZZ_p_Algebra> GFp;
typedef Field<zz_p_Algebra> gfp;
typedef Field<GF2_Algebra>  GFp2;

int main(int argv, char* argc[]) {
	double cputime;

	cout << "Using " << gfp::Infrastructure::name << endl << endl;

	cputime = -NTL::GetTime();
	const gfp* K = &(gfp::createField(3,10));
	cputime += NTL::GetTime();
	cout << *K << " in " << cputime << endl;
	cout << "Time spent building the irreducible polynomial : "
		<< gfp::TIME.BUILDIRRED << endl << endl;
	
	for (int i = 1 ; i <= 4 ; i++) {
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << *K << " in " << cputime << endl;
	}
	
	cout << endl << "Time spent building the cyclotomic polynomial : "
		<< gfp::TIME.CYCLOTOMIC << endl;


	cout << endl << endl;


	cout << "Using " << GFp2::Infrastructure::name << endl << endl;

	cputime = -NTL::GetTime();
	const GFp2* L = &(GFp2::createField(2,99));
	cputime += NTL::GetTime();
	cout << *L << " in " << cputime << endl;
	cout << "Time spent building the irreducible polynomial : "
		<< GFp2::TIME.BUILDIRRED << endl << endl;
	
	for (int i = 1 ; i <= 4 ; i++) {
		cputime = -NTL::GetTime();
		L = &(L->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << *L << " in " << cputime << endl;
	}
	
	cout << endl << "Time spent building the cyclotomic polynomial : "
		<< GFp2::TIME.CYCLOTOMIC << endl;

}
