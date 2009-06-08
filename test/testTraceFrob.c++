#include "faast.hpp"
#include <cstdlib>

using namespace std;
using namespace FAAST;

typedef Field<ZZ_p_Algebra> GFp;
typedef Field<zz_p_Algebra> gfp;
typedef Field<GF2_Algebra>  GFp2;
typedef FieldElement<ZZ_p_Algebra> GFp_E;
typedef FieldElement<zz_p_Algebra> gfp_E;
typedef FieldElement<GF2_Algebra>  GFp2_E;

int main(int argv, char* argc[]) {
	double cputime;

	gfp::Infrastructure::BigInt p;
	long d, l, t;
	cin >> p; cin >> d; cin >> l; cin >> t;

	cout << "Using " << gfp::Infrastructure::name << endl << endl;
	cputime = -NTL::GetTime();
	const gfp* K = &(gfp::createField(p,d));
	cputime += NTL::GetTime();
	cout << *K << " in " << cputime << endl;
	cout << "Time spent building the irreducible polynomial : "
		<< gfp::TIME.BUILDIRRED << endl << endl;

	for (int i = 1 ; i <= l ; i++) {
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << *K << " in " << cputime << endl;

		for (int j = 0 ; j < i+t ; j++) {
			gfp_E a = K->random(), b;
			long n = (j>=i) ?
				d + RandomBnd(K->d - d) : d*power_long(p, j);
			double frobtime, pseudotime, naivetime;

			frobtime = -NTL::GetTime();
			b = a.frobenius(n);
			frobtime += NTL::GetTime();
			cout << "Fast " << n << "-th Frobenius in "
				<< frobtime << endl;

			pseudotime = -NTL::GetTime();
			b = a.pseudotrace(n);
			pseudotime += NTL::GetTime();
			cout << "Fast " << n << "-th Pseudotrace in "
				<< pseudotime << endl;

			cout << "Pseudotraces precomputed in " <<
				gfp::TIME.PSEUDOTRACES << endl;

			naivetime = -NTL::GetTime();
			for (long i = 0 ; i < 10 ; i++)
				a.self_frobenius();
			naivetime += NTL::GetTime();
			double average = naivetime / 10;
			cout << "Naive Frobenius in " << average
				<< ", projected execution time " <<
				average * n << endl;

			cout << "Projected thresholds :" << endl
				<< "\tfrobenius " << floor(frobtime / average) << endl
				<< "\tpseudotrace " << floor(pseudotime / average) << endl;
		}
		cout << endl;
	}
}
