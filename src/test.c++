#include "Types.hpp"
#include "Field.hpp"
#include <cstdlib>

using namespace std;
using namespace AS;

typedef Field<GF2_Algebra>  gfp;
typedef FieldElement<GF2_Algebra>  gfp_E;
typedef FieldPolynomial<GF2_Algebra>  gfp_X;

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

	cout << "\tCreate\tMinPol\tInterp\tEval" << endl;
	for (int i = 1 ; i <= l ; i++) {
		cout << i << "\t";
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		gfp_E a, b, c;
		a = K->random();
		b = K->random();

		vector<gfp_X> minpols;
		cputime = -NTL::GetTime();
		a.minimalPolynomials(K->baseField(),minpols);
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		cputime = -NTL::GetTime();
		gfp_X poly = a.affineMinimalPolynomial(K->baseField(), b, minpols);
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		cputime = -NTL::GetTime();
		c = a.evaluate(poly, minpols);
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		if (c != b) {
			cout << "ERROR 1 : Results don't match" << endl;
			cout << a << endl << b << endl << c << endl;
			cout << poly << endl;
		}

		poly >>= *K;
		c = poly.evaluate(a);
		if (c != b) {
			cout << "ERROR 2 : Results don't match" << endl;
			cout << a << endl << b << endl << c << endl;
			cout << poly << endl;
		}

	}
}
