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
		
		// X mod Phi(X), the (2p-1)th root of unity
		GFpE omega; GFpX omegaX;
		SetX(omegaX); conv(omega, omegaX);
		// prepare the values X^i for i in [0..2p-2]
		vector<GFpE> omegas;
		omegas.reserve(2*p-1);
		omegas[0] = 1;
		for (long i = 1; i <= 2*p-2 ; i++)
			omegas[i] = omegas[i-1]*omega;

		// Qstar = prod_i Q(omega^i Y)
		GFpEX Qstar; conv(Qstar, Q);
		GFpEX Qtmp;
		// WARNING : oughta do better by a subproduct tree
		// approach
		for (long i = 1 ; i <= 2*p-2 ; i++) {
			long c = 0;
			for (long j = 0 ; j < deg(Q) ; j++) {
				SetCoeff(Qtmp, j,
					coeff(Q, j) * omegas[c]);
				c += i; c %= 2*p - 1;
			}
			Qstar *= Qtmp; 
		}
		
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
		
		// result = qstar(X^p - X)
		GFpX xpminusx;
		SetCoeff(xpminusx, p, 1);
		SetCoeff(xpminusx, 1, -1);
		compose<T>(res, qstar, xpminusx, p);
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
				alpha = new FieldElement<T>(one());
			} else {
				// alpha = x0
				alpha = new FieldElement<T>(*primitive);
				GFpX Q0 = GFpE::modulus().val();
				// if needed, make the trace of the primitive element
				// different from 0
				if (primitive->trace() == 0) {
					// extend modulo X^p - X - x0 - 1
					po = true;
					if (d % p == 0) throw
						NotSupportedException("I don't know how to build a primitive tower over this base field.");
					// X-1
					GFpX xminus1; SetX(xminus1); SetCoeff(xminus1,0,-1);
					// Q_0* = Q_0(X-1)
					compose<T>(Q0, Q0, xminus1, p);
					// alpha = x0 + 1
					*alpha += one();
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
			alpha = new FieldElement<T>(*primitive);
		}
		// generic case, see paper
		else {
			// extend modulo X^p - X - x1^(2p-1)
			po = false; tpmo = true;
			const GFpX& Q0 = GFpE::modulus().val();
			// switch modulus to compute modulo the cyclotomic
			// polynomial
			if ( !(baseField().Phi.get()) ) {
				GFpX phi;
				cyclotomic<T>(phi, 2*long(p)-1, p);
				GFpE::init(phi);
				Context* ctxt = new Context();
				ctxt->P.save();
				baseField().Phi.reset(ctxt);
			} else baseField().Phi->P.restore();
			// apply Cantor's algorithm to compute the
			// minimal polynomial
			cantor89<T>(Q, Q0, p);
			// alpha = x1^(2p-1)
			alpha = new FieldElement<T>(*primitive);
			*alpha ^= long(2)*p - 1;
		}
		
		// prepare the new context
		GFpE::init(Q);
		Context ctxt = baseField().context;
		ctxt.P.save();
		
		// primitive element and generator
		GFpX priX; SetX(priX);
		GFpE pri; conv(pri, priX);
		
		overfield = new Field<T>(this, ctxt, pri, po, tpmo,
								 p, long(p)*d, height+1, alpha);
		return *overfield;
	}
	
	/* Build the splitting field of the polynomial
	 * 			X^p - X - alpha
	 * over this field. This may or may not be an extension of
	 * degree p depending if the polynomial is irreducible.
	 */
//	const Field<T>& ArtinSchreierExtension(const FieldElement<T>& alpha) const throw (CharacteristicTooLargeExeption, NotSupportedException);


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
