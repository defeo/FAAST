#include "faast.hpp"
#include <cstdlib>

using namespace std;
using namespace FAAST;

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
	const gfp* K = &(gfp::createField(2,1));
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
			gfp_E a = K->random(), b;
			vector<gfp_E> down;

			cputime = -GetTime();
			pushDown(a, down);
			cputime += GetTime();
			cout << "Push-down computed in " << cputime << endl;

			cputime = -GetTime();
			liftUp(down, b);
			cputime += GetTime();
			cout << "Lift-up computed in " << cputime << endl;
			cout << "Time spent in Lift-up precomputation : " <<
				gfp::TIME.LIFTUP << endl;
			cout << "Time spent in Lift-up Transposed Multiplication : " <<
				gfp::TIME.LU_TRANSMUL << endl;
			cout << "Time spent in Lift-up Mod* : " <<
				gfp::TIME.LU_TRANSMOD << endl;
			cout << "Time spent in Lift-up Push-down-rec* : " <<
				gfp::TIME.LU_TRANSPUSHDOWN << endl;
			cout << "Time spent in Lift-up step 4 : " <<
				gfp::TIME.LU_STEP4 << endl;
			cout << "Time spent in Lift-up step 5 : " <<
				gfp::TIME.LU_STEP5 << endl;

			if (a != b) {
				cout << "ERROR : Results don't match" << endl;
				cout << a << endl << b << endl;
				vector<gfp_E>::iterator it;
				for (it = down.begin() ; it != down.end() ; it++)
					cout << *it << " ";
				cout << endl;
			}
			cout << endl;
		}
	}

	cout << endl << "Time spent building the cyclotomic polynomial : "
		<< gfp::TIME.CYCLOTOMIC << endl << endl;
}
