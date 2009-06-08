#include "faast.hpp"
#include <cstdlib>

using namespace std;
using namespace FAAST;

typedef Field<zz_p_Algebra>  gfp;
typedef FieldElement<zz_p_Algebra>  gfp_E;
typedef FieldPolynomial<zz_p_Algebra>  gfp_X;

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

		gfp_X poly;
		bool nopol = false;
		cputime = -NTL::GetTime();
		try {
			poly = a.affineMinimalPolynomial(K->baseField(), b, minpols);
		} catch (NoSuchPolynomialException e) {
			nopol = true;
			cout << "*";
		}
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		if (!nopol) {
			cputime = -NTL::GetTime();
			c = a.evaluate(poly, minpols);
			cputime += NTL::GetTime();
			cout << cputime << endl;

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
		} else cout << "*\t" << endl;
	}
}
