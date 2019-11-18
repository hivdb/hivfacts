package edu.stanford.hivdb.hivfacts;

import java.util.LinkedHashSet;
import java.util.Set;

import org.apache.commons.lang3.builder.HashCodeBuilder;

/**
 * Helper class mainly used to build mutation search index.
 *
 * Many mutation-related indices use gene and position as their index
 * key. This class instantiates hashable and comparable objects using
 * value gene and pos.
 */
public class GenePosition<U extends Gene<?, U, W>, W extends GenePosition<U, W>> implements WithGene<U>, Comparable<W> {

	public U gene;
	public Integer position;
	
	public static <T extends Strain<U>, U extends Gene<T, U, W>, W extends GenePosition<U, W>> Set<W> getGenePositionsBetween(
		final W start, final W end
	) {
		U startGene = start.getGene();
		U endGene = end.getGene();
		if (startGene.getStrain() != endGene.getStrain()) {
			throw new IllegalArgumentException(
				"Virus strain of `start` and `end` positions must be the same."
			);
		}
		T strain = startGene.getStrain();
		Set<W> genePositions = new LinkedHashSet<>();
		for (U gene : strain.getGenes()) {
			int startPos = 1;
			int endPos = gene.getLength();
			if (gene.compareTo(startGene) == 0) {
				startPos = start.getPosition();
			}
			if (gene.compareTo(endGene) == 0) {
				endPos = end.getPosition();
			}
			if (gene.compareTo(startGene) >= 0 &&
				gene.compareTo(endGene) <= 0) {
				genePositions.addAll(
					gene.getGenePositionsBetween(startPos, endPos));
			}
		}
		return genePositions;

	}

	public GenePosition(final U gene, final int position) {
		this.gene = gene;
		this.position = position;
	}

	public GenePosition(final U gene, final Integer position) {
		this.gene = gene;
		this.position = position;
	}

	public Integer getPosition() {
		return position;
	}

	public U getGene() {
		return gene;
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(63261, 362788935)
			.append(gene).append(position).toHashCode();
	}

	@Override
	public int compareTo(W o) {
		if (o == null) throw new NullPointerException("Null is incomprable.");
		int cmp = getGene().compareTo(o.getGene());
		if (cmp == 0) {
			cmp = getPosition().compareTo(o.getPosition());
		}
		return cmp;
	}

	@Override
	public String toString() {
		return String.format("%s:%d", getGene(), getPosition());
	}
}
