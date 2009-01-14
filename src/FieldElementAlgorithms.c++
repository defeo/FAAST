namespace AS {
/****************** Arithmetics ******************/
//	void self_frobenius() throw();
//	void self_frobenius(const long) throw();
//	void self_trace(const Field<T> F) throw(NotASubFieldException);
//	void self_pseudotrace(const long n) throw();
	
/****************** Level embedding ******************/
	/* Push the element e down along the stem and store
	 * the result in v.
	 * 
	 * throw : NoSubFieldException if e belongs to GF(p)
	 */
	template <class T>
	void pushDown(FieldElement<T>& e, vector<FieldElement<T> >& v)
	throw(NoSubFieldException) {
		
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
//	template <class T>
//	void liftUp(vector<FieldElement<T> >& v, FieldElement<T>& e)
//	throw(NotInSameFieldException, NoOverFieldException);
}
