package edu.stanford.hivdb.hivfacts;

public enum HIVStrain implements Strain<HIVGene> {
	HIV1("HIV-1", "hiv1b"),
	HIV2A("HIV-2 Group A", "hiv2a"),
	HIV2B("HIV-2 Group B", "hiv2b");

	private final String displayText;
	private final String nucaminoProfile;

	private HIVStrain(String displayText, String nucaminoProfile) {
		this.displayText = displayText;
		this.nucaminoProfile = nucaminoProfile;
	}

	public String getDisplayText() {
		return displayText;
	}

	public String getNucaminoProfile() {
		return nucaminoProfile;
	}

	public HIVGene[] getGenes() {
		return HIVGene.values(this);
	}
}
