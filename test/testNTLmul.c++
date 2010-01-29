#include <faast.hpp>
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
	if (cin.peek() != EOF) {
	  cin >> p; cin >> d; cin >> l;
	} else {
	  p = 2; d = 1; l = 4;
	}

	cout << "Using " << gfp::Infrastructure::name << endl << endl;
	const gfp* K = &(gfp::createField(p,d));
	cout << *K << " in " << cputime << endl << endl;
	cout << "\tFAAST\tNTL" << endl;
	for (int i = 1 ; i <= l ; i++) {
		cout << i << "\t";
		/** Construction **/
		K = &(K->ArtinSchreierExtension());

		gfp_E a = K->random(), b = K->random();

		/** FAAST Multiplication **/
		cputime = -GetTime();
		for (int i = 0 ; i < 10 ; i++)
			a*b;
		cputime += GetTime();
		cout << cputime/10 << "\t";

		/** NTL Multiplication **/
		gfp::Infrastructure::GFpE A; a.toInfrastructure(A);
		gfp::Infrastructure::GFpE B; a.toInfrastructure(B);
		gfp::Infrastructure::GFpE C;
		cputime = -GetTime();
		for (int i = 0 ; i < 10 ; i++)
			NTL::mul(C,A,B);
		cputime += GetTime();
		cout << cputime/10 << endl;
	}
}
