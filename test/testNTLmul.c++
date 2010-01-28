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
	cout << *K << " in " << cputime << endl << endl;
	cout << "\tFAAST\tNTL" << endl;
	for (int i = 1 ; i <= l ; i++) {
		cout << i << "\t";
		/** Construction **/
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << cputime << "\t";

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
		cputime = -GetTime();
		for (int i = 0 ; i < 10 ; i++)
			NTL::mul(A,B);
		cputime += GetTime();
		cout << cputime/10 << endl;
	}
	totaltime += GetTime();
}
