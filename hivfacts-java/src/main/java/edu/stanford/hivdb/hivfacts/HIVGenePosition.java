package edu.stanford.hivdb.hivfacts;

import org.apache.commons.lang3.builder.EqualsBuilder;

/**
 * Helper class mainly used to build mutation search index.
 *
 * Many mutation-related indices use gene and position as their index
 * key. This class instantiates hashable and comparable objects using
 * value gene and pos.
 */
public class HIVGenePosition extends GenePosition<HIVGene, HIVGenePosition> {

	private HIVGenePosition(String[] strainGenePos) {
		super(HIVGene.valueOf(strainGenePos[0]), Integer.parseInt(strainGenePos[1]));
	}
	
	public HIVGenePosition(final String text) {
		this(text.split(":", 2));
	}
	
	public HIVGenePosition(HIVGene gene, int pos) {
		super(gene, pos);
	}
	
	public HIVGenePosition(HIVGene gene, Integer pos) {
		super(gene, pos);
	}

	public Integer getPolPosition() {
		// internal function, don't expose to GraphQL
		int absPos;
		HIVStrain strain = gene.getStrain();
		switch(gene.getGeneEnum()) {
			case PR:
				absPos = position;
				break;
			case RT:
				absPos = HIVGene.valueOf(strain, "PR").getLength() + position;
				break;
			default:  // case IN
				absPos = (
					HIVGene.valueOf(strain, "PR").getLength() +
					HIVGene.valueOf(strain, "RT").getLength() +
					position);
				break;
		}
		return absPos;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null) { return false; }
		if (obj == this) { return true; }
		if (obj.getClass() != getClass()) { return false; }
		HIVGenePosition gp = (HIVGenePosition) obj;
		return new EqualsBuilder()
			.append(getGene(), gp.getGene())
			.append(getPosition(), gp.getPosition())
			.isEquals();
	}

}