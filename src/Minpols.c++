namespace AS {
/****************** Minimal polynomials ******************/
	/* The minimal polynomial over the field F */
	template<class T> FieldPolynomial<T>
	FieldElement<T>::minimalPolynomial(const Field<T>& F)
	const throw(NotASubFieldException) {
		if (*(F.stem) == parent_field.stem) {
			FieldPolynomial minPol = -(*this);
			minPol.setCoeff(1);
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
		Field<T>& G = parent_field->stem;
		if (!F.isSubField(G)) throw NotASubFieldException();
		
		res.resize(F.);
	}
}