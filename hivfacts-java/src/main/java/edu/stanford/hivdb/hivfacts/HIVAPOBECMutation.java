package edu.stanford.hivdb.hivfacts;

public class HIVAPOBECMutation {
	protected String gene;
	protected Integer position;
	protected Character aa;

	protected HIVAPOBECMutation(
			String gene, int position, char aa) {
		this.gene = gene;
		this.position = position;
		this.aa = aa;
	}
	
	public GenePosition getGenePosition() {
		return new GenePosition(gene, position);
	}
	
	public String getGene() { return gene; }
	public Integer getPosition() { return position; }
	public Character getAA() { return aa; }
	
	@Override
	public String toString() {
		return String.format("%s:%d%s", gene, position, aa);
	}
}