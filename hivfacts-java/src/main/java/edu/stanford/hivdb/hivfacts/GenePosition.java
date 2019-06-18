package edu.stanford.hivdb.hivfacts;

import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

public class GenePosition implements Comparable<GenePosition> {
	protected final String gene;
	protected final Integer position;
	
	protected GenePosition(final String gene, final int pos) {
		this.gene = gene;
		this.position = pos;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj == null) { return false; }
		if (obj == this) { return true; }
		if (obj.getClass() != getClass()) { return false; }
		GenePosition gp = (GenePosition) obj;
		return new EqualsBuilder()
			.append(gene, gp.gene)
			.append(position, gp.position)
			.isEquals();
	}
	
	@Override
	public int compareTo(GenePosition o) {
		if (o == null) throw new NullPointerException("Null is incomprable.");
		int cmp = gene.compareTo(o.gene);
		if (cmp == 0) {
			cmp = Integer.valueOf(position).compareTo(o.position);
		}
		return cmp;
	}
	
	@Override
	public String toString() {
		return String.format("%s:%d", gene, position);
	}
	
	@Override
	public int hashCode() {
		return new HashCodeBuilder(4663259, 1290575637)
			.append(gene).append(position).toHashCode();
	}

}
