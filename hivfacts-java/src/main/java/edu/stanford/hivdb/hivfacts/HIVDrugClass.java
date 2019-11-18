package edu.stanford.hivdb.hivfacts;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import edu.stanford.hivdb.hivfacts.HIVGeneEnum;

public enum HIVDrugClass implements DrugClass<HIVDrug> {
	NRTI(HIVGeneEnum.RT, "Nucleoside Reverse Transcriptase Inhibitor"),
	NNRTI(HIVGeneEnum.RT, "Non-nucleoside Reverse Transcriptase Inhibitor"),
	PI(HIVGeneEnum.PR, "Protease Inhibitor"),
	INSTI(HIVGeneEnum.IN, "Integrase Strand Transfer Inhibitor");

	private final HIVGeneEnum gene;
	private final String fullName;

	private HIVDrugClass(final HIVGeneEnum gene, final String fullName) {
		this.gene = gene;
		this.fullName = fullName;
	}

	public List<HIVDrug> getAllDrugs() {
		return Stream.of(HIVDrug.values())
			   .filter(d -> d.getDrugClass() == this)
			   .collect(Collectors.toList());
	}

	public HIVGeneEnum gene() {
		return this.gene;
	}

	public String getFullName() {
		return this.fullName;
	}

	public static HIVDrugClass getSynonym(String synonym) {
		Map<String, HIVDrugClass> drugSynonyms = new HashMap<>();
		drugSynonyms.put("INI", HIVDrugClass.INSTI);
		if (!drugSynonyms.containsKey(synonym)) {
			return HIVDrugClass.valueOf(synonym);
		} else {
			return drugSynonyms.get(synonym);
		}
	}
}
