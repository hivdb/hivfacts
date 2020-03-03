package edu.stanford.hivdb.hivfacts;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import edu.stanford.hivdb.comments.ConditionalComments;
import edu.stanford.hivdb.drugresistance.algorithm.DrugResistanceAlgorithm;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.genotypes.Genotype;
import edu.stanford.hivdb.genotypes.GenotypeReference;
import edu.stanford.hivdb.genotypes.Genotyper;
import edu.stanford.hivdb.mutations.AminoAcidPercents;
import edu.stanford.hivdb.mutations.CodonPercents;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationPrevalence;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.mutations.MutationTypePair;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

public class HIV implements Virus<HIV> {

	private static final String VIRUS_NAME = "HIV";
	private static final String MAIN_STRAIN = "HIV1";
	private static final String STRAINS_RESPATH = "strains_hiv1.json";
	private static final String GENES_RESPATH = "genes_hiv1.json";
	private static final String DRUG_CLASSES_RESPATH = "drug-classes_hiv1.json";
	private static final String DRUGS_RESPATH = "drugs.json";
	private static final String DRMS_RESPATH = "drms_hiv1.json";
	private static final String SDRMS_RESPATH = "sdrms_hiv1.json";
	private static final String TSMS_RESPATH = "tsms_hiv1.json";
	private static final String APOBECS_RESPATH = "apobecs/apobecs.json";
	private static final String APOBEC_DRMS_RESPATH = "apobecs/apobec_drms.json";
	private static final String AAPCNTS_RESPATH = "aapcnt/rx-%s_subtype-%s.json";
	private static final String CODONPCNTS_RESPATH = "codonpcnt/rx-%s_subtype-%s.json";
	private static final String MUTTYPES_RESPATH = "mutation-types.json";
	private static final String MUTTYPE_PAIRS_RESPATH = "mutation-type-pairs_hiv1.json";
	private static final String MAIN_SUBTYPES_RESPATH = "main-subtypes_hiv1.json";
	private static final String GENOTYPE_REFERENCES_RESPATH = "genotypes/genotype-references_hiv1.json";
	private static final String GENOTYPES_RESPATH = "genotypes/genotypes_hiv1.json";
	private static final String ALGORITHMS_INDEXPATH = "algorithms/versions.json";
	private static final String ALGORITHMS_RESPATH = "algorithms/%s_%s.xml";
	private static final String CONDCOMMENTS_RESPATH = "conditional-comments_hiv1.json";

	static {
		Virus.registerInstance(new HIV());
	}
		
	public static HIV getInstance() {
		return Virus.getInstance(HIV.class);
	}
	
	private final HIVDataLoader<HIV> dl;

	private HIV() {
		registerSequenceValidator(new HIVDefaultSequenceValidator());
		registerMutationsValidator(new HIVDefaultMutationsValidator());
		registerSequenceReadsValidator(new HIVDefaultSequenceReadsValidator());
		this.dl = new HIVDataLoader<>(
			this,
			VIRUS_NAME,
			MAIN_STRAIN,
			STRAINS_RESPATH,
			GENES_RESPATH,
			DRUG_CLASSES_RESPATH,
			DRUGS_RESPATH,
			DRMS_RESPATH,
			SDRMS_RESPATH,
			TSMS_RESPATH,
			APOBECS_RESPATH,
			APOBEC_DRMS_RESPATH,
			AAPCNTS_RESPATH,
			CODONPCNTS_RESPATH,
			MUTTYPES_RESPATH,
			MUTTYPE_PAIRS_RESPATH,
			MAIN_SUBTYPES_RESPATH,
			GENOTYPE_REFERENCES_RESPATH,
			GENOTYPES_RESPATH,
			ALGORITHMS_INDEXPATH,
			ALGORITHMS_RESPATH,
			CONDCOMMENTS_RESPATH
		);
	}

	@Override
	public String getName() {
		return dl.getName();
	}
	
	@Override
	public Strain<HIV> getMainStrain() {
		return dl.getMainStrain();
	}
	
	@Override
	public Collection<Strain<HIV>> getStrains() {
		return dl.getStrains();
	}
	
	@Override
	public Strain<HIV> getStrain(String name) {
		return dl.getStrain(name);
	}

	@Override
	public Collection<Gene<HIV>> getGenes(Strain<HIV> strain) {
		return dl.getGenes(strain);
	}
	
	@Override
	public Gene<HIV> getGene(String name) {
		return dl.getGene(name);
	}
	
	@Override
	public Collection<DrugClass<HIV>> getDrugClasses() {
		return dl.getDrugClasses();
	}
	
	@Override
	public Map<String, DrugClass<HIV>> getDrugClassSynonymMap() {
		return dl.getDrugClassSynonymMap();
	}
	
	@Override
	public DrugClass<HIV> getDrugClass(String name) {
		return dl.getDrugClass(name);
	}
	
	@Override
	public Collection<Drug<HIV>> getDrugs() {
		return dl.getDrugs();
	}
	
	@Override
	public Map<String, Drug<HIV>> getDrugSynonymMap() {
		return dl.getDrugSynonymMap();
	}

	@Override
	public Collection<DrugResistanceAlgorithm<HIV>> getDrugResistAlgorithms() {
		return dl.getDrugResistAlgorithms();
	}

	@Override
	public Collection<DrugResistanceAlgorithm<HIV>> getDrugResistAlgorithms(Collection<String> algorithmNames) {
		return dl.getDrugResistAlgorithms(algorithmNames);
	}
	
	
	@Override
	public DrugResistanceAlgorithm<HIV> getDrugResistAlgorithm(String name) {
		return dl.getDrugResistAlgorithm(name);
	}

	@Override
	public DrugResistanceAlgorithm<HIV> getDrugResistAlgorithm(String family, String version) {
		return dl.getDrugResistAlgorithm(family, version);
	}
	
	@Override
	public Gene<HIV> extractMutationGene(String mutText) {
		return dl.extractMutationGene(mutText);
	}

	@Override
	public Mutation<HIV> parseMutationString(Gene<HIV> defaultGene, String mutText) {
		return dl.parseMutationString(defaultGene, mutText);
	}

	@Override
	public Mutation<HIV> parseMutationString(String mutText) {
		return dl.parseMutationString(mutText);
	}
	
	@Override
	public MutationSet<HIV> newMutationSet(String formattedMuts) {
		return dl.newMutationSet(formattedMuts);
	}

	@Override
	public MutationSet<HIV> newMutationSet(Collection<String> formattedMuts) {
		return dl.newMutationSet(formattedMuts); 
	}

	@Override
	public MutationSet<HIV> newMutationSet(Gene<HIV> defaultGene, String formattedMuts) {
		return dl.newMutationSet(defaultGene, formattedMuts);
	}

	@Override
	public MutationSet<HIV>	newMutationSet(Gene<HIV> defaultGene, Collection<String> formattedMuts) {
		return dl.newMutationSet(defaultGene, formattedMuts);
	}

	@Override
	public Map<DrugClass<HIV>, MutationSet<HIV>> getDrugResistMutations() {
		return dl.getDrugResistMutations();
	}
	
	@Override
	public Map<DrugClass<HIV>, MutationSet<HIV>> getSurveilDrugResistMutations() {
		return dl.getSurveilDrugResistMutations();
	}

	@Override
	public Map<DrugClass<HIV>, MutationSet<HIV>> getRxSelectedMutations() {
		return dl.getRxSelectedMutations();
	}
	
	@Override
	public MutationSet<HIV> getApobecMutations() {
		return dl.getApobecMutations();
	}

	@Override
	public MutationSet<HIV> getApobecDRMs() {
		return dl.getApobecDRMs();
	}

	@Override
	public Collection<MutationType<HIV>> getMutationTypes() {
		return dl.getMutationTypes();
	}
	
	@Override
	public MutationType<HIV> getMutationType(String mutTypeText) {
		return dl.getMutationType(mutTypeText);
	}

	@Override
	public Collection<MutationTypePair<HIV>> getMutationTypePairs() {
		return dl.getMutationTypePairs();
	}
	
	@Override
	public AminoAcidPercents<HIV> getAminoAcidPercents(Strain<HIV> strain, String treatment, String subtype) {
		return dl.getAminoAcidPercents(strain, treatment, subtype);
	}

	@Override
	public CodonPercents<HIV> getCodonPercents(Strain<HIV> strain, String treatment, String subtype) {
		return dl.getCodonPercents(strain, treatment, subtype);
	}

	@Override
	public List<MutationPrevalence<HIV>> getMutationPrevalence(GenePosition<HIV> genePos) {
		return dl.getMutationPrevalence(genePos);
	}
	
	@Override
	public ConditionalComments<HIV> getConditionalComments() {
		return dl.getConditionalComments();
	}
	
	@Override
	public List<String> getMainSubtypes(Strain<HIV> strain) {
		return dl.getMainSubtypes(strain);
	}
	
	@Override
	public Map<Gene<HIV>, Map<String, Integer[]>> getNumPatientsForAAPercents(Strain<HIV> strain) {
		return dl.getNumPatientsForAAPercents(strain);
	}

	@Override
	public Collection<Genotype<HIV>> getGenotypes() {
		return dl.getGenotypes();
	}
	
	@Override
	public Genotype<HIV> getGenotype(String name) {
		return dl.getGenotype(name);
	}

	@Override
	public Genotype<HIV> getGenotypeUnknown() {
		return dl.getGenotypeUnknown();
	}

	@Override
	public List<GenotypeReference<HIV>> getGenotypeReferences() {
		return dl.getGenotypeReferences();
	}
	
	@Override
	public Genotyper<HIV> getGenotyper() {
		return dl.getGenotyper();
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// HIV instance is a singleton
		return false;
	}

}