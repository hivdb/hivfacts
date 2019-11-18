package edu.stanford.hivdb.hivfacts;

import java.util.List;

public interface DrugClass<T extends Drug<?>> {

	public List<T> getAllDrugs();

	public String getFullName();

}
