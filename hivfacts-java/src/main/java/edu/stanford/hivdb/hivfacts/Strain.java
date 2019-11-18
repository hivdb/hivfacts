package edu.stanford.hivdb.hivfacts;

public interface Strain<U extends Gene<?, U, ?>> {

	public String getDisplayText();

	public String getNucaminoProfile();

	public U[] getGenes();

}
