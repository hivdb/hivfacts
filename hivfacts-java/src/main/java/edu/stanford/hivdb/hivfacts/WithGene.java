package edu.stanford.hivdb.hivfacts;

public interface WithGene<U extends Gene<?, U, ?>> {
	
	public U getGene();

}
