#include "utilities.hpp"

namespace AS {
	template <class T> typename T::MatGFp artinMatrix
	(const typename T::BigInt& p, const long line, const typename T::GFpXModulus& P) {
		typedef typename T::MatGFp      MatGFp;
		typedef typename T::GFpX        GFpX;

		long d = deg(P);
		// we just want an invertinle d-1 minor
		MatGFp appl; appl.SetDims(d-1, d-1);
		// X^p mod P
		GFpX Xp = PowerXMod(p, P);
		// column is the column of X^p-X at each iteration
		GFpX column(0,1);
		// build the minor
		for (long i = 1 ; i < d ; i++) {
			MulMod(column, column, Xp, P);
			for (long j = 0 ; j < line ; j++) {
				appl[j][i-1] = coeff(column, j);
			}
			for (long j = line+1 ; j < d ; j++) {
				appl[j-1][i-1] = coeff(column, j);
			}
			if (i > line) appl[i-1][i-1]--;
			else if (i < line) appl[i][i-1]--;
		}
		//invert
		return inv(appl);
	}

/****************** Access to precomputed values ******************/
	template <class T> const FieldElement<T>&
	Field<T>::getPseudotrace(const long j) const {
#ifdef AS_DEBUG
		if (j < 0 || j >= height)
			throw ASException("Bad input to getPseudotrace.");
#endif

		long size = pseudotraces.size();
		if (size > j) return pseudotraces[j];
		if (this != stem) return stem->getPseudotrace(j);

#ifdef AS_TIMINGS
		TIME.PSEUDOTRACES = -GetTime();
#endif
		pseudotraces.resize(j+1);
		if (size == 0) {
			pseudotraces[0] = *alpha;
			pseudotraces[0].SmallPTrace(baseField().d);
			size++;
		}
		for (long i = size ; i <= j ; i++) {
			pseudotraces[i] = pseudotraces[i-1];
			FieldElement<T> t = pseudotraces[i];
			for (BigInt h = 1 ; h < p ; h++) {
				t.BigFrob(i-1);
				pseudotraces[i] += t;
			}
		}
#ifdef AS_TIMINGS
		TIME.PSEUDOTRACES += GetTime();
#endif

		return pseudotraces[j];
	}
	
	template <class T> const FieldElement<T>& Field<T>::getLiftup() const {
		if (this != stem) return stem->getLiftup();
		
		if ( !(liftuphelper.get()) ) {
			switchContext();
#ifdef AS_TIMINGS
			TIME.LIFTUP = -GetTime();
#endif
			GFpX diffQ; diff(diffQ, GFpE::modulus());
			FieldElement<T>* helper = new FieldElement<T>();
			helper->base = false;
			conv(helper->repExt, diffQ);
			helper->parent_field = this;
			
			helper->self_inv();
			liftuphelper.reset(helper);
#ifdef AS_TIMINGS
			TIME.LIFTUP += GetTime();
#endif
		}
		return *liftuphelper;
	}
	
	template <class T> const typename
	T::MatGFp& Field<T>::getArtinMatrix() const {
		if (this != stem) return stem->getArtinMatrix();

		if (artinLine == -1) {
			// We pick a redundant line : it corresponds
			// to a power of x of trace different from 0.
			// The residue formula tells us that Tr(x^dep) != 0
			switchContext();
			GFpXModulus P = GFpE::modulus();
#ifdef AS_TIMINGS
			TIME.ARTINMATRIX = -GetTime();
#endif
			artinLine = d - 1 - deg(diff(P));
			artin = artinMatrix<T>(p,artinLine,P);
#ifdef AS_TIMINGS
			TIME.ARTINMATRIX += GetTime();
#endif
		}
		return artin;
	}

	template <class T> const typename
	T::Context& Field<T>::getCyclotomic() const {
		if (this != stem) return stem->getCyclotomic();

		switchContext();
		if ( !(Phi.get()) ) {
			GFpX phi;
#ifdef AS_TIMINGS
			TIME.CYCLOTOMIC = -GetTime();
#endif
			cyclotomic<T>(phi, 2*long(p)-1, p);
#ifdef AS_TIMINGS
			TIME.CYCLOTOMIC += GetTime();
#endif
			GFpE::init(phi);
			Context* ctxt = new Context();
			ctxt->P.save();
			Phi.reset(ctxt);
		}
		return *Phi;
	}

}
