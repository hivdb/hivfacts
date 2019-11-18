package edu.stanford.hivdb.hivfacts;

public class HIVAminoAcidPercent implements WithGene<HIVGene> {
	final private String strain = "HIV1";
	final private String gene;
	final public Integer position;
	final public Character aa;
	final public Double percent;
	final public Integer count;
	final public Integer total;
	final public String reason;
	final public Boolean isUnusual;

	protected HIVAminoAcidPercent (
			String gene, int position, char aa, double percent,
			int count, int total, String reason, boolean isUnusual) {
		this.gene = gene;
		this.position = position;
		this.aa = aa;
		this.percent = percent;
		this.count = count;
		this.total = total;
		this.reason = reason;
		this.isUnusual = isUnusual;
	}

	@Override
	public HIVGene getGene() {
		return HIVGene.valueOf(strain, gene);
	}

	public HIVGenePosition getGenePosition() {
		return new HIVGenePosition(getGene(), position);
	}

}
