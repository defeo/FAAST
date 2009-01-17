#include "utilities.hpp"
#include "Tmul.hpp"

NTL_OPEN_NNS
void ShiftAdd(GF2X& c, const GF2X& a, long n);
void ShiftAdd(zz_pX& c, const zz_pX& a, long n);
void ShiftAdd(ZZ_pX& c, const ZZ_pX& a, long n);
void ComputeTraceVec(const GF2XModulus& F);
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

	/* The routine Transposed mul from Section 4
	 * (step 2 of Lift-up)
	 * Computes the transposed multiplication of the linear form
	 *     (0, ..., 0, -form)
	 * by the element
	 *     (W[0], ..., W[p-1])
	 * modulo Q
	 */
	template <class T> void TransposedMul(
	vector<typename T::GFpX>& W, const typename T::GFpXModulus& Q,
	const typename T::GFpX& form, const typename T::BigInt& p) {
		typedef typename T::GFpX            GFpX;
		typedef typename T::GFpXMultilplier GFpXMultiplier;
		typedef typename T::BigInt          BigInt;

		vector<GFpXMultiplier> Trans; Trans.resize(p);
		
		build(Trans[p-1], W[0] + W[long(p)-1], Q);
		for (BigInt i = 1; i < p ; i++)
			build(Trans[p-i-1], W[i], Q);
		
		GFpX formtmp = -form;
		for (BigInt i = 0 ; i < p ; i++)
			TransMulMod(W[i], formtmp, Trans[i], Q);
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
		if (!e.parent_field->subfield)
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
			// the result lies in GF(p)[x0+1].
			// This brings the elements back to GF(p)[x0]
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
			e.parent_field->subfield->switchContext();
			
			// convert the result of push-down-rec to elements
			// of the subfield
			v.resize(p);
			bool base = e.parent_field->subfield->d == 1;
			for (BigInt i = 0 ; i < p ; i++) {
				v[i].base = base;
				// this automatically reduces modulo
				// the defining polynomial if needed
				if (base) {
					v[i].repBase = coeff(W[i], 0);
					v[i].repExt = 0;
				} else {
					v[i].repBase = 0;
					conv(v[i].repExt, W[i]);
				}
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
		typedef typename T::GFpX        GFpX;
		typedef typename T::GFpE        GFpE;
		typedef typename T::GFpXModulus GFpXModulus;
		typedef typename T::BigInt      BigInt;

		// if v is empty, return the generic 0
		if (v.size() == 0) {
			e = FieldElement<T>();
			return;
		}
		// check all elements of v belong to the same field
		// and save this field into parent
		const Field<T>* parent = v[0].parent_field;
		for (long i = 0 ; i < v.size() ; i++) {
			if (!parent) parent = v[i].parent_field;
			else if (v[i].parent_field != parent)
				throw NotInSameFieldException();
		}
		// standard checks
		if (!parent) throw NoOverFieldException();
		if (parent->stem != parent) throw NoOverFieldException();
		if (!parent->overfield) throw NoOverFieldException();
		
		parent->switchContext();
		
		// if this is the prime field of a base field
		// simply merge the coefficients
		if (parent->overfiled->height == 0) {
			GFpX eX;
			for (long i = min(parent->overfield->d, v.size()) - 1 ; i >=0 ; i--) {
				if (!v[i].isZero())
					SetCoeff(eX, i, v[i].repBase());
			}
			parent->overfield->switchContext();
			e.base = false;
			e.repBase = 0;
			conv(e.repExt, eX);
			e.parent_field = parent->overfield;
		}
		// the real lift-up algorithm from Section 4
		else {
			bool base = parent->d == 1;
			BigInt p = parent->p;

			// take the elements out of v
			vector<GFpX> W; W.resize(p);
			for (BigInt i = 0 ; i < p ; i++) {
				if (!v[i].isZero()) {
					if (base) W[i] = v[i].repBase;
					else W[i] = rep(v[i].repExt());
				}
			}
			
			// The input lies in GF(p)[x0].
			// If this extension was built modulo
			//   X^p - X - x0 - 1
			// this brings the elements into GF(p)[x0+1]
			if (parent->overfield->plusone) {
				GFpX xminusone;
				SetCoeff(xminusone, 1); SetCoeff(xminusone, 0, -1);
				for (BigInt i = 0 ; i < p ; i++)
					compose<T>(W[i], W[i], xminusone, p);
			}

			// get the trace form
			const GFpXModulus& Q = GFpE::modulus();
			if (Q.tracevec.length() == 0)
				ComputeTraceVec(Q);
			GFpX trace; conv(trace, Q.tracevec);
			// TransposedMul (step 2 of lift-up)
			TransposedMul<T>(W, Q, trace, p);

			// if this extension was built modulo
			//   X^p - X - xi^(2p-1)
			// apply mod* and evaluate*
			// (steps 2 and 3 of push-down*)
			if (parent->overfield->twopminusone) {
				// mod*
				
				
				// evaluate*
				for (BigInt i = 0 ; i < p ; i++)
					contract<T>(W[i], W[i], 2*long(p) - 1);
			}
			
			// step 4 of push-down*
			//TransPushDownRec();
			
			//
		}
	}
	
}
