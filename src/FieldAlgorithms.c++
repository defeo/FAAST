#include "utilities.hpp"

/* This file contains algorithms from Sections 3 and 6
 * of the paper
 */

namespace AS {
	
	/* Cantor's algorithm to compute the minimal polynomial
	 * of the primitive generator of the Artin-Schreier
	 * extension.
	 * 
	 * It assumes the modulus has already been properly set
	 * to the (2p-1)th cyclotomic polynomial.
	 */
	template <class T> void cantor89(typename T::GFpX& res,
	const typename T::GFpX& Q, const long p) {
		typedef typename T::GFpX  GFpX;
		typedef typename T::GFpE  GFpE;
		typedef typename T::GFpEX GFpEX;
		
#ifdef AS_TIMINGS
		Field<T>::TIME.C89_PRE = -GetTime();
#endif
		// X mod Phi(X), the (2p-1)th root of unity
		GFpE omega; GFpX omegaX;
		SetX(omegaX); conv(omega, omegaX);
		// prepare the values X^i for i in [0..2p-2]
		GFpE omegas[2*p-1];
		omegas[0] = 1;
		for (long i = 1; i <= 2*p-2 ; i++)
			omegas[i] = omegas[i-1]*omega;
#ifdef AS_TIMINGS
		Field<T>::TIME.C89_PRE += GetTime();
		Field<T>::TIME.C89_Qstar = -GetTime();
#endif

		// Qstar = prod_i Q(omega^i Y)
		GFpEX Qstar; conv(Qstar, Q);
		GFpEX Qtmp;
		// WARNING : oughta do better by a subproduct tree
		// approach
		for (long i = 1 ; i <= 2*p-2 ; i++) {
			long c = 0;
			SetCoeff(Qtmp, deg(Q));  // hack to speed up GF2
			for (long j = 0 ; j <= deg(Q) ; j++) {
				SetCoeff(Qtmp, j,
					coeff(Q, j) * omegas[c]);
				c += i; c %= 2*p - 1;
			}
			Qstar *= Qtmp; 
		}
#ifdef AS_TIMINGS
		Field<T>::TIME.C89_Qstar += GetTime();
		Field<T>::TIME.C89_qstar = -GetTime();
#endif
		// qstar( X^(2p-1) ) = Qstar
		GFpX qstar;
		long c = 0;
		for (long i = 0 ; i <= deg(Qstar) ; i += 2*p-1) {
			SetCoeff(qstar, c, coeff(rep(coeff(Qstar, i)), 0));
			c++;
#ifdef AS_DEBUG
			if (deg(rep(coeff(Qstar, i))) > 0) throw
				ASException("Error in Cantor89");
			for (long j = 1 ; j < 2*p-1 ; j++) {
				if (coeff(Qstar, i+j) != 0) throw
				ASException("Error in Cantor89");
			}
#endif
		}
#ifdef AS_TIMINGS
		Field<T>::TIME.C89_qstar += GetTime();
		Field<T>::TIME.C89_compose = -GetTime();
#endif
		
		// result = qstar(X^p - X)
		GFpX xpminusx;
		SetCoeff(xpminusx, p, 1);
		SetCoeff(xpminusx, 1, -1);
		compose<T>(res, qstar, xpminusx, p);
#ifdef AS_TIMINGS
		Field<T>::TIME.C89_compose += GetTime();
#endif
	}


/****************** Field Extensions ******************/
	/* Build a default extension of degree p over this field. 
	 */
	template <class T> const Field<T>& Field<T>::ArtinSchreierExtension()
	const throw (CharacteristicTooLargeException, NotSupportedException) {
#ifdef AS_DEBUG
		if (!stem) throw ASException("Error : Stem is NULL.");
#endif
		// if the extension already exists, return it
		// WARNING : this behavior is not correct when
		// this field is the prime field of some other field
		if (stem->overfield) return *(stem->overfield);

		// Build the extension (Section 3)

		// test if the characteristic stays in one word
		if (p != long(p)) throw CharacteristicTooLargeException();
		
#ifdef AS_TIMINGS
		TIME.BUILDSTEM = -GetTime();
#endif
		switchContext();
		GFpX Q; bool po, tpmo;
		FieldElement<T>* alpha;

		// Compute the minimal polynomial Q of the new level
		 
		// first floor, Q_1 = Q_0(X^p - X)
		if (height == 0) {
			// extend modulo X^p - X - x0
			tpmo = false; po = false;
			// if this field is GF(p)
			if (d == 1) {
				// the minimal polynomials is X^p - X - 1
				SetCoeff(Q, p, 1);
				SetCoeff(Q, 1, -1);
				SetCoeff(Q, 0, -1);
				// alpha = 1
				alpha = new FieldElement<T>(stem->one());
			} else {
				// alpha = x0
				alpha = new FieldElement<T>(*(stem->primitive));
				GFpX Q0 = GFpE::modulus().val();
				// if needed, make the trace of the primitive element
				// different from 0
				if (stem->primitive->trace() == 0) {
					// extend modulo X^p - X - x0 - 1
					po = true;
					if (d % p == 0) throw
						NotSupportedException("I don't know how to build a primitive tower over this base field.");
					// X-1
					GFpX xminus1; SetX(xminus1); SetCoeff(xminus1,0,-1);
					// Q_0* = Q_0(X-1)
					compose<T>(Q0, Q0, xminus1, p);
					// alpha = x0 + 1
					*alpha += stem->one();
				}
				// X^p - X
				GFpX xpminusx;
				SetCoeff(xpminusx, p, 1);
				SetCoeff(xpminusx, 1, -1);
				// Q_1 = Q_0
				compose<T>(Q, Q0, xpminusx, p);
			}
		}
		// second floor in characteristic 2,
		// Q_2 = Q_1(X^p - X)
		else if (height == 1 && p == long(2)) {
			// extend modulo X^p - X - x1
			po = false; tpmo = false;
			const GFpX& Q0 = GFpE::modulus().val();
			// X^p - X
			GFpX xpminusx;
			SetCoeff(xpminusx, p, 1);
			SetCoeff(xpminusx, 1, -1);
			// Q_1 = Q_0
			compose<T>(Q, Q0, xpminusx, p);
			// alpha = x1
			alpha = new FieldElement<T>(*(stem->primitive));
		}
		// generic case, see paper
		else {
			// extend modulo X^p - X - x1^(2p-1)
			po = false; tpmo = true;
			const GFpX& Q0 = GFpE::modulus().val();
			// switch modulus to compute modulo the cyclotomic
			// polynomial
			primeField().getCyclotomic().P.restore();
			// apply Cantor's algorithm to compute the
			// minimal polynomial
#ifdef AS_TIMINGS
				TIME.CANTOR89 = -GetTime();
#endif
			cantor89<T>(Q, Q0, p);
#ifdef AS_TIMINGS
				TIME.CANTOR89 += GetTime();
#endif
			// alpha = x1^(2p-1)
			alpha = new FieldElement<T>(*(stem->primitive));
			*alpha ^= long(2)*p - 1;
		}
		
#if AS_DEBUG >= 3
		if (!IterIrredTest(Q))
			throw ASException("The defining polynomial of the extension is not irreducible.");
#endif
		// prepare the new context
		GFpE::init(Q);
		Context ctxt = primeField().context;
		ctxt.P.save();
		
		// primitive element and generator
		GFpX priX; SetX(priX);
		GFpE pri; conv(pri, priX);
		
		// who generated this extension ?
		const Field<T>* vsub = (stem == this)? NULL : this;
		
		stem->overfield = new Field<T>(stem, ctxt, pri, po, tpmo, p,
										long(p)*d, height+1, alpha, vsub);
#ifdef AS_TIMINGS
		TIME.BUILDSTEM += GetTime();
#endif

		return *(stem->overfield);
	}
	
	/* Build the splitting field of the polynomial
	 * 			X^p - X - alpha
	 * over this field. This may or may not be an extension of
	 * degree p depending if the polynomial is irreducible.
	 */
	template <class T> const Field<T>&
	Field<T>::ArtinSchreierExtension(const FieldElement<T>& alpha)
	const throw (CharacteristicTooLargeException,
	NotSupportedException, IllegalCoercionException) {
		// check that alpha belongs to here
		if ( !isOverFieldOf(*(alpha.parent_field)) )
			throw IllegalCoercionException();
		if (alpha.isZero()) return *this;
		
		FieldElement<T> root;
		FieldElement<T>* aleph = new FieldElement<T>(alpha);
		const Field<T> *vsub, *st;
		// if X^p - X - alpha generates an extension of degree p
		if (alpha.parent_field->stem == stem && alpha.trace() != 0) {
			// first build the stem if needed
			const Field<T>& up = ArtinSchreierExtension();
			// find a root in the overfield
			root = up.Couveignes2000(alpha);
			vsub = this; st = up.stem;
		}
		// if X^p - X - alpha generates an extension of degree 1
		else {
			// find a root in here
			root = Couveignes2000(alpha);
			vsub = stem; st = stem;
		}
		// move alpha in here
		*aleph >>= *this;
		return *(new Field<T>(st, root, aleph, vsub));
	}


/****************** Level embedding ******************/
	/* Push the element e down to this field and store
	 * the result in v.
	 * 
	 * throw : IllegalCoercionException if the field e belongs to
	 *         is not the immediate overfield of this.
	 */
//	void pushDown(FieldElement<T>& e, vector<FieldElement<T> >& v) const
//		 throw(IllegalCoercionException);

	/* Lift the elements in v up to this field and store the result in e.
	 * If v is too short, it is filled with zeros. If v is too
	 * long, the unnecessary elements are ignored.
	 * 
	 * throw : NotInSameFieldException if the elements of v do not
	 *         belong all to the same field.
	 * throw : IllegalCoercionException if the field e belongs to
	 *         is not the immediate subfield of this.
	 */
//	void liftUp(vector<FieldElement<T> >& v, FieldElement<T>& e) const
//		throw(NotInSameFieldException, IllegalCoercionException);

}
