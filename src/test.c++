#include "Types.hpp"
#include "Field.hpp"
#include <cstdlib>

using namespace std;
using namespace AS;

typedef Field<ZZ_p_Algebra> GFp;
typedef Field<zz_p_Algebra> gfp;
typedef Field<GF2_Algebra>  GFp2;

int main(int argv, char* argc[]) {
	const GFp& L = GFp::createField(to_ZZ(3),5);
	const gfp& k = gfp::createField(2,9);
	const GFp2& K = GFp2::createField(2,9);
	const GFp2& base = K.baseField();

	cout << K.characteristic() << " " << K.degree() << endl 
		<< k.characteristic() << " " << k.degree() << endl
		<< L << ", " << base << endl
		<< L.cardinality() << " " << K.cardinality() << " "
		<< base.cardinality() << " " << k.cardinality() << endl;
	cout << K.isOverFieldOf(base) << K.isIsomorphic(K) << (K==K) << endl;
	
	double cputime = -NTL::GetTime();
	const GFp2& KK = K.ArtinSchreierExtension();
	cputime += NTL::GetTime();
	cout << KK << " in " << cputime << endl;
	
	cputime = -NTL::GetTime();
	const GFp2& KKK = KK.ArtinSchreierExtension();
	cputime += NTL::GetTime();
	cout << KKK << " in " << cputime << endl;
	
	cputime = -NTL::GetTime();
	const GFp2& KKKK = KKK.ArtinSchreierExtension();
	cputime += NTL::GetTime();
	cout << KKKK << " in " << cputime << endl;
	
	cputime = -NTL::GetTime();
	const GFp2& KKKKK = KKKK.ArtinSchreierExtension();
	cputime += NTL::GetTime();
	cout << KKKKK << " in " << cputime << endl;
}
