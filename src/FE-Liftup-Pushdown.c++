#include "utilities.hpp"
#include "Tmul.hpp"

NTL_OPEN_NNS
void ShiftAdd(GF2X& c, const GF2X& a, long n);
void ShiftAdd(zz_pX& c, const zz_pX& a, long n);
void ShiftAdd(ZZ_pX& c, const ZZ_pX& a, long n);
void ComputeTraceVec(const GF2XModulus& F);
NTL_CLOSE_NNS

namespace AS {
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
		
		GFpX Lead; RightShift(Lead, W[0], 1);
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
				pushDownRec<T>(V, i, min(i+splitdegree-1, end), Wtmp, p);
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
				for (BigInt j = 0 ; j < p ; j++)
					Wtmp[j] = trunc(W[j], splitdegree/p);
				TransMulMod<T>(W, k-1, p);
				TransPushDownRec<T>(Wtmp, V, i, min(i+splitdegree-1, end), p);
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
		typedef typename T::GFpX           GFpX;
		typedef typename T::GFpXMultiplier GFpXMultiplier;
		typedef typename T::BigInt         BigInt;

		vector<GFpXMultiplier> Trans; Trans.resize(p);
		
		build(Trans[long(p)-1], W[0] + W[long(p)-1], Q);
		for (BigInt i = 1; i < p ; i++)
			build(Trans[long(p)-long(i)-1], W[i], Q);
		
		GFpX formtmp = -form;
		for (BigInt i = 0 ; i < p ; i++)
			TransMulMod(W[i], formtmp, Trans[i], Q);
	}

	template <class T> void TransMod(typename T::GFpX& W,
	const typename T::GFpXModulus& Q, const typename T::BigInt& p) {
		typedef typename T::GFpX           GFpX;
		typedef typename T::GFpXMultiplier GFpXMultiplier;
		typedef typename T::BigInt         BigInt;
		
		long d = deg(Q);
		// Xn = X^d mod Q
		GFpX Xn = -Q; SetCoeff(Xn, d, 0);
		GFpXMultiplier Trans; build(Trans, Xn, Q);

		GFpX tmp1 = W;
		GFpX tmp2;
		long shift = 0;
		for (BigInt i = 1 ; i < 2*long(p)-1 ; i++) {
			TransMulMod(tmp2, tmp1, Trans, Q);
			tmp1 = tmp2;
			shift += d;
			ShiftAdd(W, tmp1, shift);
		}
	}
	
	
	
	/* Push the element e down along the stem and store
	 * the result in v.
	 * 
	 * throw : NoSubFieldException if e belongs to GF(p)
	 */
	template <class T>
	void pushDown(const FieldElement<T>& e, vector<FieldElement<T> >& v)
	throw(NoSubFieldException) {
		typedef typename T::GFp    GFp;
		typedef typename T::GFpX   GFpX;
		typedef typename T::BigInt BigInt;

		if (!e.parent_field) {
			v.clear();
			return;
		}
		const Field<T>* parent = e.parent_field->stem;
		if (!parent->subfield)
			throw NoSubFieldException();
		
		parent->switchContext();
		
		// if the subfield is prime
		// simply return the list of coefficients
		if (parent->subfield->d == 1) {
			v.resize(deg(rep(e.repExt)) + 1);
			const Field<T>* base = &(parent->primeField());
			const GFpX& eX = rep(e.repExt);
			for (long i = 0 ; i <= deg(eX) ; i++) {
				v[i].base = true;
				v[i].repBase = coeff(eX, i);
				v[i].repExt = 0;
				v[i].parent_field = base;
			}
		}
		// if the element is a scalar
		else if (deg(rep(e.repExt)) <= 0) {
			v.resize(1);
			v[0].base = false;
			v[0].repBase = 0;
			GFp ec = coeff(rep(e.repExt), 0);
			parent->subfield->switchContext();
			v[0].repExt = ec;
			v[0].parent_field = parent->subfield;
		}
		// the real push-down algorithm from Section 4
		else {
			vector<GFpX> W;
			BigInt p = parent->p;
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
			if (parent->plusone) {
				GFpX xplusone;
				SetCoeff(xplusone, 1); SetCoeff(xplusone, 0);
				for (BigInt i = 0 ; i < p ; i++)
					compose<T>(W[i], W[i], xplusone, p);
			}
			// if this extension was built modulo
			//   X^p - X - xi^(2p-1)
			if (parent->twopminusone) {
				for (BigInt i = 0 ; i < p ; i++)
					expand<T>(W[i], W[i], 2*long(p) - 1);
			}
			
			// prepare to work in the subfield
			parent->subfield->switchContext();
			
			// convert the result of push-down-rec to elements
			// of the subfield
			v.resize(p);
			bool base = parent->subfield->d == 1;
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
				v[i].parent_field = parent->subfield;
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
	void liftUp(const vector<FieldElement<T> >& v, FieldElement<T>& e)
	throw(NotInSameFieldException, NoOverFieldException) {
		typedef typename T::GFp    GFp;
		typedef typename T::GFpX        GFpX;
		typedef typename T::GFpE        GFpE;
		typedef typename T::GFpXModulus GFpXModulus;
		typedef typename T::BigInt      BigInt;

		// if v is empty, return the generic 0
		typename vector<FieldElement<T> >::const_iterator it = v.begin();
		if (it == v.end()) {
			e = FieldElement<T>();
			return;
		}
		// check all elements of v belong to the same field
		// and save this field into parent.
		// In parallel, check if v represents a scalar
		const Field<T>* parent = it->parent_field;
		bool scalar = it->base || deg(rep(it->repExt)) <= 0;
		for (it++ ; it != v.end() ; it++) {
			if (!parent) parent = it->parent_field;
			else if (it->parent_field && it->parent_field != parent)
				throw NotInSameFieldException();
			scalar = scalar && *it == 0;
		}
		// standard checks
		if (!parent) {
			e = FieldElement<T>();
			return;
		}
		parent = parent->stem;
		if (!parent->overfield) throw NoOverFieldException();
		
		parent->switchContext();
		
		// if this is a prime field
		// simply merge the coefficients
		if (parent->d == 1) {
			GFpX eX;
			for (long i = min(parent->overfield->d, long(v.size())) - 1 ; i >=0 ; i--) {
				if (!v[i].isZero())
					SetCoeff(eX, i, v[i].repBase);
			}
			parent->overfield->switchContext();
			e.base = false;
			e.repBase = 0;
			conv(e.repExt, eX);
			e.parent_field = parent->overfield;
		}
		// if the element is a scalar
		else if (scalar) {
			e.base = false;
			e.repBase = 0;
			GFp ec = coeff(rep(v[0].repExt), 0);
			parent->overfield->switchContext();
			e.repExt = ec;
			e.parent_field = parent->overfield;
		}
		// the real lift-up algorithm from Section 4
		else {
			bool base = parent->d == 1;
			BigInt p = parent->p;
			GFpXModulus Q = GFpE::modulus();

			// take the elements out of v
			vector<GFpX> W; W.resize(p);
			for (BigInt i = 0 ; i < p ; i++) {
				if (i < long(v.size()) && !v[i].isZero()) {
					if (base) W[i] = v[i].repBase;
					else W[i] = rep(v[i].repExt);
				}
			}
			// The input lies in GF(p)[x0].
			// If this extension was built modulo
			//   X^p - X - x0 - 1
			// this brings the elements into GF(p)[x0+1]
			if (parent->overfield->plusone) {
#ifdef AS_TIMINGS
				Field<T>::TIME.LU_PLUSONE = -GetTime();
#endif
				GFpX xminusone;
				SetCoeff(xminusone, 1); SetCoeff(xminusone, 0, -1);
				for (BigInt i = 0 ; i < p ; i++)
					compose<T>(W[i], W[i], xminusone, p);
				GFpX tmp;
				compose<T>(tmp, Q.val(), xminusone, p);
				build(Q, tmp);
#ifdef AS_TIMINGS
				Field<T>::TIME.LU_PLUSONE += GetTime();
#endif
			}

			// get the trace form
			if (Q.tracevec.length() == 0) {
#ifdef AS_TIMINGS
				Field<T>::TIME.TRACEVEC = -GetTime();
#endif
				ComputeTraceVec(Q);
#ifdef AS_TIMINGS
				Field<T>::TIME.TRACEVEC += GetTime();
#endif
			}
			GFpX trace; conv(trace, Q.tracevec);
			// TransposedMul (step 2 of lift-up)
#ifdef AS_TIMINGS
			Field<T>::TIME.LU_TRANSMUL = -GetTime();
#endif
			TransposedMul<T>(W, Q, trace, p);
#ifdef AS_TIMINGS
			Field<T>::TIME.LU_TRANSMUL += GetTime();
#endif

			// if this extension was built modulo
			//   X^p - X - xi^(2p-1)
			// apply mod* and evaluate*
			// (steps 2 and 3 of push-down*)
			if (parent->overfield->twopminusone) {
				// mod*
#ifdef AS_TIMINGS
				Field<T>::TIME.LU_TRANSMOD = -GetTime();
#endif
				for (BigInt i = 0 ; i < p ; i++)
					TransMod<T>(W[i], Q, p);
#ifdef AS_TIMINGS
				Field<T>::TIME.LU_TRANSMOD += GetTime();
#endif
				// evaluate*
#ifdef AS_TIMINGS
				Field<T>::TIME.LU_TRANSEVAL = -GetTime();
#endif
				for (BigInt i = 0 ; i < p ; i++)
					contract<T>(W[i], W[i], 2*long(p) - 1);
#ifdef AS_TIMINGS
				Field<T>::TIME.LU_TRANSEVAL += GetTime();
#endif
			}
			
			// step 4 of push-down*
			GFpX V;
#ifdef AS_TIMINGS
			Field<T>::TIME.LU_TRANSPUSHDOWN = -GetTime();
#endif
			TransPushDownRec<T>(W, V, 0, parent->overfield->d - 1, p);
#ifdef AS_TIMINGS
			Field<T>::TIME.LU_TRANSPUSHDOWN += GetTime();
#endif
			
			// now get ready to work in the overfield
			parent->overfield->switchContext();
			
			// step 4 of lift-up
			const GFpXModulus& QQ = GFpE::modulus();
#ifdef AS_TIMINGS
			Field<T>::TIME.LU_STEP4 = -GetTime();
#endif
			GFpX revQ; reverse(revQ, QQ);
			MulTrunc(V, V, revQ, deg(QQ));
#ifdef AS_TIMINGS
			Field<T>::TIME.LU_STEP4 += GetTime();
#endif
			
			// step 5 of lift-up
#ifdef AS_TIMINGS
			Field<T>::TIME.LU_STEP5 = -GetTime();
#endif
			reverse(V, V, deg(QQ) - 1);
			e.base = false;
			e.repBase = 0;
			conv(e.repExt, V);
			e.parent_field = parent->overfield;
#ifdef AS_TIMINGS
			Field<T>::TIME.LU_STEP5 += GetTime();
#endif
			const FieldElement<T>& invDiffQQ =
				parent->overfield->getLiftup(); 
#ifdef AS_TIMINGS
			Field<T>::TIME.LU_STEP5 -= GetTime();
#endif
			e *= invDiffQQ;
#ifdef AS_TIMINGS
			Field<T>::TIME.LU_STEP5 += GetTime();
#endif
		}
	}
	
}
