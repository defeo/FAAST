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

	cout << "Using " << gfp::Infrastructure::name << endl << endl;
	cputime = -NTL::GetTime();
	const gfp* K = &(gfp::createField(7,1));
	cputime += NTL::GetTime();
	cout << *K << " in " << cputime << endl;
	cout << "Time spent building the irreducible polynomial : "
		<< gfp::TIME.BUILDIRRED << endl << endl;
	
	for (int i = 1 ; i <= 5 ; i++) {
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << *K << " in " << cputime << endl;

		for (int i = 1 ; i <= 3 ; i++) {
			gfp_E a = K->random(), b, c, tmp;
			long n = RandomBnd(K->d);
			cputime = -NTL::GetTime();
			b = a.frobenius(n);
			cputime += NTL::GetTime();
			cout << "Fast " << n << "-th Frobenius in "
				<< cputime << endl;
			c = a;
			cputime = -NTL::GetTime();
			for (long i = 0 ; i < n ; i++)
				c.self_frobenius();
			cputime += NTL::GetTime();
			cout << "Naive " << n << "-th Frobenius in "
				<< cputime << endl;
			if (b != c) {
				cout << "ERROR : results don't match" << endl;
				cout << a << endl << b << endl << c << endl;
			}
			cout << "Pseudotraces precomputed in " <<
				gfp::TIME.PSEUDOTRACES << endl; 
			cputime = -NTL::GetTime();
			b = a.pseudotrace(n);
			cputime += NTL::GetTime();
			cout << "Fast " << n << "-th Pseudotrace in "
				<< cputime << endl;
			c = a;
			cputime = -NTL::GetTime();
			tmp = c;
			for (long i = 1 ; i < n ; i++) {
				tmp.self_frobenius();
				c += tmp;
			}
			cputime += NTL::GetTime();
			cout << "Naive " << n << "-th Pseudotrace in "
				<< cputime << endl;
			if (b != c) {
				cout << "ERROR : results don't match" << endl;
				cout << a << endl << b << endl << c << endl;
			}
		}
		cout << endl;
	}
}
