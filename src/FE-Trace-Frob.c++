namespace AS {
/****************** Arithmetics ******************/
	template <class T> void FieldElement<T>::self_frobenius()
	throw() {
		if (!parent_field) return;
		if (base) return;
		
		parent_field->switchContext();
		
		power(repExt, repExt, parent_field->p);
	}
	
	template <class T> void
	FieldElement<T>::self_frobenius(const long) throw() {
		if (!parent_field) return;
		if (base) return;
		
		parent_field->switchContext();
		
		//TODO
	}
	
//	void self_trace(const Field<T> F) throw(NotASubFieldException);
//	void self_pseudotrace(const long n) throw();
	
}
