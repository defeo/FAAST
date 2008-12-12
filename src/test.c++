#include "Types.hpp"
#include "Field.hpp"

using namespace std;
using namespace AS;

typedef Field<ZZ_p_Algebra> GFp;
typedef Field<zz_p_Algebra> gfp;
typedef Field<GF2_Algebra>  GFp2;

int main(int argv, char* argc[]) {
	const GFp& K = GFp::createField(to_ZZ(5),4);
	const gfp& k = gfp::createField(2,5);
	const GFp2& L = GFp2::createField(2,3);
	const GFp& base = K.baseField();

	cout << K.characteristic() << " " << K.degree() << endl 
		<< k.characteristic() << " " << k.degree() << endl
		<< L << ", " << base << endl
		<< L.cardinality() << " " << K.cardinality() << " "
		<< base.cardinality() << " " << k.cardinality() << endl;
	cout << K.isOverFieldOf(base) << K.isIsomorphic(K) << (K==K) << endl;
}
