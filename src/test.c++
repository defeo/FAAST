#include "Types.h"

using namespace std;
using namespace AS;

template <class T> class boh {
public:
	void boh1 () {
		typename T::GFp p;
	}
};

int main(int argv, char* argc[]) {
//	boh<ZZ_p_Algebra> alg;
//	ZZ_pX p = alg.boh1();
	boh<ZZ_p_Algebra>();
}
