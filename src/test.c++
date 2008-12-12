#include "Types.hpp"
#include "Field.hpp"

using namespace std;
using namespace AS;

typedef Field<ZZ_p_Algebra> GFp;
typedef Field<GF2_Algebra>  GFp2;

int main(int argv, char* argc[]) {
	GFp& K = GFp::createField(5,4);
	GFp2& L = GFp2::createField(2,3);
	cout << K.characteristic() << " " << K.degree() << " " 
		<< L.characteristic() << " " << L.degree() << " "
		<< L.cardinality() << endl;
}
