namespace AS {
/****************** Minimal polynomials ******************/
	/* The minimal polynomial over the field F */
	template<class T> FieldPolynomial<T>
	FieldElement<T>::minimalPolynomial(const Field<T>& F)
	const throw(NotASubFieldException) {
		if (*(F.stem) == parent_field.stem) {
			FieldPolynomial<T> minPol = -(*this);
			minPol.setCoeff(1);
			minPol >>= *(F.stem);
			return minPol;
		} else if (*(F.stem) == parent_field.primeField()) {
			parent_field.switchContext();
			GFpX minpol;
			MinPolyMod(minpol, rep(repExt), GFpE::modulus());
			return F.stem->fromInfrastructure(minpol);
		} else {
			vector<FieldPolynomial<T> > minpols;
			minpols = minimalPolynomials(F);
			return minpols[0];
		}
	}

	/* All the minimal polynomials up to the field F.
	 *		res = minimalPolynomials(F);
	 * res[0] contains minimalPolynomials(F), res[1] contains
	 * minimalPolynomials(F.overfield) and so on.
	 */
	template <class T> void FieldElement<T>::minimalPolynomials(
	const Field<T>& F, vector<FieldPolynomial<T> > res)
	const throw(NotASubFieldException) {
		Field<T>* G = parent_field->stem;
		if (!F.isSubField(*G)) throw NotASubFieldException();
		
		long levels = G->height - F.height + 1;
		long prime = 0;
		if (F.stem->overfield->isBaseField() && !G->isPrimeField())
			prime = 1;
		res.resize(levels+prime);
		// minimal polynomial is  (X-this)
		levels--;
		res[levels + prime] = -(*this);
		res[levels + prime].setCoeff(1);
		res[levels + prime] >>= *G;
		// go down as long as this is a member of the subfield
		while (isCoercible(*(G->subfield))) {
			res[levels + prime - 1] =
				res[levels + prime] >> *(G->subfield);
			levels--;
			G = G->subfield;
		}
		// go down by galois conjugation
		while (!G->isBaseField() && *G != *(F.stem)) {
			FieldPolynomial<T> tmp = res[levels + prime];
			res[levels + prime - 1] = res[levels + prime];
			for (BigInt i = 1 ; i < G->p ; i++) {
				tmp.frobenius(power_long(G->p, G->height - 1)*G->d);
				res[levels + prime -1] *= tmp;
			}
			res[levels + prime - 1] >>= *(G->subfield);
			G = G->subfield;
			levels--;
		}
		// last level
		if (prime > 0) {
			res[0] = minimalPolynomial(F);
		}
	}
}
