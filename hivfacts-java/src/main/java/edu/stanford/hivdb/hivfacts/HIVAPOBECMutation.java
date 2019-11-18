package edu.stanford.hivdb.hivfacts;

public class HIVAPOBECMutation implements WithGene<HIVGene> {
	final private String strain = "HIV1";
	final private String gene;
	final protected Integer position;
	final protected Character aa;

	protected HIVAPOBECMutation(
			String gene, int position, char aa) {
		this.gene = gene;
		this.position = position;
		this.aa = aa;
	}

	public HIVGenePosition getGenePosition() {
		return new HIVGenePosition(getGene(), position);
	}

	@Override
	public HIVGene getGene() { return HIVGene.valueOf(strain, gene); }
	public Integer getPosition() { return position; }
	public Character getAA() { return aa; }

	@Override
	public String toString() {
		return String.format("%s%s:%d%s", strain, gene, position, aa);
	}
}
