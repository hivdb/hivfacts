package edu.stanford.hivdb.hivfacts;

public interface Drug<T extends DrugClass<?>> {
	
	public T getDrugClass();

	public String getFullName();

	public String getDisplayAbbr();

}
