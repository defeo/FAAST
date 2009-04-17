#include "Types.hpp"
#include "Field.hpp"
#include <cstdlib>

using namespace std;
using namespace AS;

typedef Field<GF2_Algebra> gfp;
typedef FieldElement<GF2_Algebra> gfp_E;

int main(int argv, char* argc[]) {
	double cputime, totaltime;

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
	
	cout << "\tCreate\tPushDow\tLiftUp\tPreLift\tMul" << endl;
	totaltime = -GetTime();
	for (int i = 1 ; i <= l ; i++) {
		cout << i << "\t";
		/** Construction **/
		cputime = -NTL::GetTime();
		K = &(K->ArtinSchreierExtension());
		cputime += NTL::GetTime();
		cout << cputime << "\t";

		gfp_E a = K->random(), b;
		vector<gfp_E> down;

		/** Push-down / Lift-up **/
		cputime = -GetTime();
		pushDown(a, down);
		cputime += GetTime();
		cout << cputime << "\t";

		cputime = -GetTime();
		liftUp(down, b);
		cputime += GetTime();
		cout << cputime << "\t";
		cout << gfp::TIME.LIFTUP << "\t";

		if (a != b) {
			cout << "ERROR : Results don't match" << endl;
			cout << a << endl << b << endl;
			vector<gfp_E>::iterator it;
			for (it = down.begin() ; it != down.end() ; it++)
				cout << *it << " ";
			cout << endl;
		}

		/** Multiplication **/
		cputime = -GetTime();
		a*b;
		cputime += GetTime();
		cout << cputime << endl;

		/** Iterated frobenius **/
/*		long n = RandomBnd(K->d - K->subField().d) + K->subField().d;
		cputime = -GetTime();
		a.self_frobenius(n);
		cputime += GetTime();
		cout << n << "-th frobenius computed in " << cputime << endl;
		cout << "Pseudotraces precomputed in " <<
			gfp::TIME.PSEUDOTRACES << endl; 
*/		
	}
	totaltime += GetTime();
		
	cout << endl << "Time spent building the cyclotomic polynomial : "
		<< gfp::TIME.CYCLOTOMIC << endl << endl;
	cout << "Total duration : " << totaltime << endl;
}
