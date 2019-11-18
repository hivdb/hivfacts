package edu.stanford.hivdb.hivfacts;

public class HIVCodonPercent implements WithGene<HIVGene> {
	final private String strain = "HIV1";
	final private String gene;
	final public Integer position;
	final public String codon;
	final public Character aa;
	final public Double percent;
	final public Integer count;
	final public Integer total;

	protected HIVCodonPercent(
			String gene, int position, String codon, char aa,
			double percent, int count, int total) {
		this.gene = gene;
		this.position = position;
		this.codon = codon;
		this.aa = aa;
		this.percent = percent;
		this.count = count;
		this.total = total;
	}

	@Override
	public HIVGene getGene() {
		return HIVGene.valueOf(strain, gene);
	}

	public HIVGenePosition getGenePosition() {
		return new HIVGenePosition(getGene(), position);
	}
}
