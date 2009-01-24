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
cout << gfp::Infrastructure::GFpE::modulus() << endl;
	cout << *K << " in " << cputime << endl;
	cout << "Time spent building the irreducible polynomial : "
		<< gfp::TIME.BUILDIRRED << endl << endl;
	
	for (int i = 1 ; i <= l ; i++) {
		gfp_E alpha = K->random();
		cout << "alpha = " << alpha << endl;
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension(alpha));
		cputime += NTL::GetTime();
cout << gfp::Infrastructure::GFpE::modulus() << endl;
		cout << *K << " in " << cputime << endl;
		cout << "Time spent building the stem " <<
			gfp::TIME.BUILDSTEM << endl;

		cout << "Tr(alpha) = " << alpha.trace() << ", root = "
			<< K->Couveignes2000(alpha) << endl << endl;
	}
	cout << "Time spent inverting the matrix" <<
		gfp::TIME.ARTINMATRIX << endl;
}
