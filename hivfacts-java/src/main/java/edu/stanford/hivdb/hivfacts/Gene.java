package edu.stanford.hivdb.hivfacts;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * enum for Gene. Will replace the previous class-wrapper-around enum implementation now called Gene
 *
 */
public abstract class Gene<T extends Strain<U>, U extends Gene<T, U, W>, W extends GenePosition<U, W>> implements Comparable<U> {

	public static final Character WILDCARD = '.';

	/**
	 * Get the drug classes associated with a gene
	 */
	public abstract List<HIVDrugClass> getDrugClasses();

	public Set<W> getGenePositionsBetween(
		int start, int end
	) {
		Set<W> genePositions = new LinkedHashSet<>();
		start = Math.max(start, 1);
		end = Math.min(end, getLength());
		for (int pos = start; pos <= end; pos ++) {
			genePositions.add(getGenePosition(pos));
		}
		return genePositions;
	}
	
	public abstract W getGenePosition(Integer pos);

	/**
	 * Get the mutation types associated with a gene
	 */
	public abstract List<MutType> getMutationTypes();

	/**
	 * Get the scored mutation types associated with a gene
	 */
	public abstract List<MutType> getScoredMutTypes();

	public abstract int getLength();

	public abstract int getFirstNA();

	public abstract int getNASize();

	/**
	 * Get the reference amino acid (AA) at a position in a gene
	 * Indexed starting from 1
	 * @param pos
	 * @return the AA at the submitted position
	 */
	public abstract String getReference(int pos);

	public abstract Character getRefChar(int pos);

	public abstract String getReference(int pos, int length);

	public abstract String getReference();

	public abstract T getStrain();

	public abstract String getNameWithStrain();

	public abstract String getName();

	@Override
	public String toString() {
		return getNameWithStrain();
	}

	@Override
	public abstract int compareTo(U o);

	@Override
	public abstract int hashCode();

	@Override
	public abstract boolean equals(Object o);
}
