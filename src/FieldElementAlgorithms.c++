#include "utilities.hpp"
#include "Tmul.hpp"

NTL_OPEN_NNS
void ShiftAdd(GF2X& c, const GF2X& a, long n);
void ShiftAdd(zz_pX& c, const zz_pX& a, long n);
void ShiftAdd(ZZ_pX& c, const ZZ_pX& a, long n);
NTL_CLOSE_NNS

namespace AS {
/****************** Arithmetics ******************/
//	void self_frobenius() throw();
//	void self_frobenius(const long) throw();
//	void self_trace(const Field<T> F) throw(NotASubFieldException);
//	void self_pseudotrace(const long n) throw();
	
/****************** Level embedding ******************/
	// The routine MulMod from Section 4
	template <class T> void MulMod(
	vector<typename T::GFpX>& W, const long n,
	const typename T::BigInt& p) {
		typedef typename T::GFpX   GFpX;
		typedef typename T::BigInt BigInt;
		
		GFpX Lead = W[long(p)-1];
		for (BigInt i = p-long(1) ; i >= long(0) ; i--) {
			GFpX tmp = W[i];
			W[i] = 0;
			long shift = 1;
			for (long j = 0 ; j < n ; j++) {
				ShiftAdd(W[i], tmp, shift);
				shift *= p;
			}
			if (i > long(0)) W[i] += W[long(i)-1];
		}
		W[1] += Lead;
		ShiftAdd(W[0], Lead, 1);
	}

	// The routine MulMod* from Section 4
	template <class T> void TransMulMod(
	vector<typename T::GFpX>& W, const long n,
	const typename T::BigInt& p) {
		typedef typename T::GFpX   GFpX;
		typedef typename T::BigInt BigInt;
		
		GFpX Lead = RightShift(W[0], 1);
		Lead += W[1];
		for (BigInt i = 0 ; i <= p-long(1) ; i++) {
			if (i > long(0)) W[long(i)-1] += W[i];
			GFpX tmp;
			long shift = 1;
			for (long j = 0 ; j < n ; j++) {
				tmp += RightShift(W[i], shift);
				shift *= p;
			}
			W[i] = tmp;
		}
		W[long(p)-1] += Lead;
	}

	// The routine Push-down-rec from Section 4
	template <class T> void pushDownRec(
	const typename T::GFpX& V, long start, long end,
	vector<typename T::GFpX>& W, const typename T::BigInt& p) {
		typedef typename T::GFpX   GFpX;
		typedef typename T::BigInt BigInt;

		long degree = end - start;
		long k = NumPits(p, degree);
		// if deg(V) >= p, cut in p slices and apply recursively
		if (k > 1) {
			vector<GFpX> Wtmp;
			W.clear(); W.resize(p);
			long splitdegree = power_long(p, k-1);
			for (long i = start + splitdegree * (degree / splitdegree) ; i >= start ; i -= splitdegree) {
				pushDownRec<T>(V, i, min(i+splitdegree-1, deg(V)), Wtmp, p);
				MulMod<T>(W, k-1, p);
				for (BigInt j = 0 ; j < p ; j++) W[j] += Wtmp[j];
			}
		}
		// if deg(V) < p, then
		//     V mod Z^p - Z - T
		// is V 
		else {
			W.resize(p);
			for (long i = start ; i <= end ; i++)
				W[i-start] = coeff(V, i);
		}
	}

	// The routine Push-down-rec* from Section 4
	template <class T> void TransPushDownRec(
	vector<typename T::GFpX>& W, typename T::GFpX& V,
	long start, long end, const typename T::BigInt& p) {
		typedef typename T::GFpX   GFpX;
		typedef typename T::BigInt BigInt;

		long degree = end - start;
		long k = NumPits(p, degree);
		// if deg(V) >= p, cut in p slices and apply recursively
		if (k > 1) {
			vector<GFpX> Wtmp; Wtmp.resize(p);
			long splitdegree = power_long(p, k-1);
			SetCoeff(V, end); // hack
			for (long i = start ; i <= end ; i += splitdegree) {
				for (BigInt j = 0 ; j < p ; j++) Wtmp[j] += W[j];
				TransMulMod<T>(W, k-1, p);
				TransPushDownRec<T>(Wtmp, V, i, min(i+splitdegree-1, deg(V)), p);
			}
		}
		else {
			for (long i = end ; i >= start ; i--)
				SetCoeff(V, i, coeff(W[i-start], 0));
		}
	}

	// The routine Push-down* from Section 4
	template <class T> void PowerProjection(
	vector<typename T::GFpX>& W, typename T::GFpX& V,
	const typename T::GFpX& Q, long d, const typename T::BigInt& p) {
		typedef typename T::GFpX   GFpX;
		typedef typename T::BigInt BigInt;
		
		for (BigInt i = 0; i < p ; i++) TransMod(W[i], Q, 2*p-1);
		
		// if this extension was built modulo
		//   X^p - X - xi^(2p-1)
		if (twopminusone) {
			for (BigInt i = 0 ; i < p ; i++)
				contract<T>(W[i], W[i], 2*long(p) - 1);
		}
		
		// if this extension was built modulo
		//   X^p - X - x0 - 1
		if (plusone) {
			GFpX xplusone;
			SetCoeff(xplusone, 1); SetCoeff(xplusone, 0);
			for (BigInt i = 0 ; i < p ; i++)
				TransCompose<T>(W[i], W[i], xplusone, p);
		}
		
		TransPushDownRec<T>(W, V, 0, d, p);
	}



	/* Push the element e down along the stem and store
	 * the result in v.
	 * 
	 * throw : NoSubFieldException if e belongs to GF(p)
	 */
	template <class T>
	void pushDown(const FieldElement<T>& e, vector<FieldElement<T> >& v)
	throw(NoSubFieldException) {
		typedef typename T::GFpX   GFpX;
		typedef typename T::BigInt BigInt;

		if (!e.parent_field)
			throw NoSubFieldException();
		if (e.parent_field->stem != e.parent_field)
			throw NoSubFieldException();
		if (e.parent_field->d == 1)
			throw NoSubFieldException();
		
		e.parent_field->switchContext();
		
		// if this is a (non-prime) base field
		// simply return the list of coefficients
		if (e.parent_field->height == 0) {
			v.resize(e.parent_field->d);
			const Field<T>* base = &(e.parent_field->baseField());
			const GFpX& eX = rep(e.repExt);
			for (long i = 0 ; i < e.parent_field->d ; i++) {
				v[i].base = true;
				v[i].repBase = coeff(eX, i);
				v[i].repExt = 0;
				v[i].parent_field = base;
			}
		}
		// the real push-down algorithm from Section 4
		else {
			vector<GFpX> W;
			BigInt p = e.parent_field->p;
			pushDownRec<T>(rep(e.repExt), 0, deg(rep(e.repExt)), W, p);
#if AS_DEBUG >= 2
			for (BigInt i = 0 ; i < p ; i++) {
				if (deg(W[i]) * p > deg(rep(e.repExt)))
					throw ASException("Problem in pushDownRec.");
			}
#endif
			// if this extension was built modulo
			//   X^p - X - x0 - 1
			// (could do better implenting this directly
			// into pushDownRec)
			if (e.parent_field->plusone) {
				GFpX xplusone;
				SetCoeff(xplusone, 1); SetCoeff(xplusone, 0);
				for (BigInt i = 0 ; i < p ; i++)
					compose<T>(W[i], W[i], xplusone, p);
			}
			// if this extension was built modulo
			//   X^p - X - xi^(2p-1)
			if (e.parent_field->twopminusone) {
				for (BigInt i = 0 ; i < p ; i++)
					expand<T>(W[i], W[i], 2*long(p) - 1);
			}
			
			// prepare to work in the subfield
#ifdef AS_DEBUG
			if (!e.parent_field->subfield)
				throw ASException("In pushDown : no subfield.");
#endif
			e.parent_field->subfield->switchContext();
			
			// convert the result of push-down-rec to elements
			// of the subfield
			v.resize(p);
			for (BigInt i = 0 ; i < p ; i++) {
				v[i].base = false;
				v[i].repBase = 0;
				// this automatically reduces modulo
				// the defining polynomial
				conv(v[i].repExt, W[i]);
				v[i].parent_field = e.parent_field->subfield;
			}
		}
	}
		
	/* Lift the elements in v up along the stem and store
	 * the result in e.
	 * If v is too short, it is filled with zeros.
	 * If v is too long, the unnecessary elements are ignored.
	 * 
	 * throw : NotInSameFieldException if the elements of v do not
	 *         belong all to the same field.
	 * throw : NoOverFieldException if there's no extension to lift
	 *         up to.
	 */
	template <class T>
	void liftUp(const vector<const FieldElement<T> >& v, FieldElement<T>& e)
	throw(NotInSameFieldException, NoOverFieldException) {
		typedef typename T::GFpX GFpX;

		
	}
}
