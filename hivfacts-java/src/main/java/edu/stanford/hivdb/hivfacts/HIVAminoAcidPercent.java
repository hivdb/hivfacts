package edu.stanford.hivdb.hivfacts;

public class HIVAminoAcidPercent {
	final public String gene;
	final public Integer position;
	final public Character aa;
	final public Double percent;
	final public Integer count;
	final public Integer total;
	final public String reason;
	final public Boolean isAPOBEC;
	final public Boolean isUnusual;

	protected HIVAminoAcidPercent(
			String gene, int position, char aa, double percent, int count,
			int total, String reason, boolean isAPOBEC, boolean isUnusual) {
		this.gene = gene;
		this.position = position;
		this.aa = aa;
		this.percent = percent;
		this.count = count;
		this.total = total;
		this.reason = reason;
		this.isAPOBEC = isAPOBEC;
		this.isUnusual = isUnusual;
	}
	
	public GenePosition getGenePosition() {
		return new GenePosition(gene, position);
	}
}
