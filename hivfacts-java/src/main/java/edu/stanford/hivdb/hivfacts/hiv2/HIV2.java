package edu.stanford.hivdb.hivfacts.hiv2;

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
import edu.stanford.hivdb.hivfacts.HIVDataLoader;
import edu.stanford.hivdb.mutations.AminoAcidPercents;
import edu.stanford.hivdb.mutations.CodonPercents;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationPrevalence;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.mutations.MutationTypePair;
import edu.stanford.hivdb.seqreads.SequenceReadsAssembler;
import edu.stanford.hivdb.sequences.AlignmentConfig;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

public class HIV2 implements Virus<HIV2> {

	private static final String VIRUS_NAME = "HIV2";
	private static final String MAIN_STRAIN = "HIV2A";
	private static final String STRAINS_RESPATH = "strains_hiv2.json";
	private static final String GENES_RESPATH = "genes_hiv2.json";
	private static final String DRUG_CLASSES_RESPATH = "drug-classes_hiv2.json";
	private static final String DRUGS_RESPATH = "drugs.json";
	private static final String DRMS_RESPATH = "drms_hiv2.json";
	private static final String SDRMS_RESPATH = "sdrms_hiv2.json";
	private static final String TSMS_RESPATH = "tsms_hiv2.json";
	private static final String APOBECS_RESPATH = "apobecs-hiv2/apobecs.json";
	private static final String APOBEC_DRMS_RESPATH = "apobecs-hiv2/apobec_drms.json";
	private static final String AAPCNTS_RESPATH = "aapcnt-hiv2/rx-%s_subtype-%s.json";
	private static final String CODONPCNTS_RESPATH = "codonpcnt-hiv2/rx-%s_subtype-%s.json";
	private static final String MUTTYPES_RESPATH = "mutation-types.json";
	private static final String MUTTYPE_PAIRS_RESPATH = "mutation-type-pairs_hiv2.json";
	private static final String MAIN_SUBTYPES_RESPATH = "main-subtypes_hiv2.json";
	private static final String GENOTYPE_REFERENCES_RESPATH = "genotypes/genotype-references_hiv2.json";
	private static final String GENOTYPES_RESPATH = "genotypes/genotypes_hiv2.json";
	private static final String ALGORITHMS_INDEXPATH = "algorithms/versions_hiv2.json";
	private static final String ALGORITHMS_RESPATH = "algorithms/%s-HIV2_%s.xml";
	private static final String CONDCOMMENTS_RESPATH = "conditional-comments_hiv2.json";
	private static final String ALIGNCONFIG_RESPATH = "alignment-config_hiv2.json";
	private static final String ASSEMBLYCONFIG_RESPATH = "assembly-config_hiv2.json";

	static {
		Virus.registerInstance(new HIV2());
	}
		
	public static HIV2 getInstance() {
		return Virus.getInstance(HIV2.class);
	}

	private final HIVDataLoader<HIV2> dl;
	
	private HIV2() {
		registerSequenceValidator(new HIV2DefaultSequenceValidator());
		registerMutationsValidator(new HIV2DefaultMutationsValidator());
		registerSequenceReadsValidator(new HIV2DefaultSequenceReadsValidator());
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
			CONDCOMMENTS_RESPATH,
			ALIGNCONFIG_RESPATH,
			ASSEMBLYCONFIG_RESPATH
		);
	}

	@Override
	public String getName() {
		return dl.getName();
	}
	
	@Override
	public Strain<HIV2> getMainStrain() {
		return dl.getMainStrain();
	}
	
	@Override
	public Collection<Strain<HIV2>> getStrains() {
		return dl.getStrains();
	}
	
	@Override
	public Strain<HIV2> getStrain(String name) {
		return dl.getStrain(name);
	}

	@Override
	public Collection<Gene<HIV2>> getGenes(Strain<HIV2> strain) {
		return dl.getGenes(strain);
	}
	
	@Override
	public Gene<HIV2> getGene(String name) {
		return dl.getGene(name);
	}
	
	@Override
	public Collection<DrugClass<HIV2>> getDrugClasses() {
		return dl.getDrugClasses();
	}
	
	@Override
	public Map<String, DrugClass<HIV2>> getDrugClassSynonymMap() {
		return dl.getDrugClassSynonymMap();
	}
	
	@Override
	public DrugClass<HIV2> getDrugClass(String name) {
		return dl.getDrugClass(name);
	}
	
	@Override
	public Collection<Drug<HIV2>> getDrugs() {
		return dl.getDrugs();
	}
	
	@Override
	public Map<String, Drug<HIV2>> getDrugSynonymMap() {
		return dl.getDrugSynonymMap();
	}

	@Override
	public Collection<DrugResistanceAlgorithm<HIV2>> getDrugResistAlgorithms() {
		return dl.getDrugResistAlgorithms();
	}

	@Override
	public Collection<DrugResistanceAlgorithm<HIV2>> getDrugResistAlgorithms(Collection<String> algorithmNames) {
		return dl.getDrugResistAlgorithms(algorithmNames);
	}
	
	
	@Override
	public DrugResistanceAlgorithm<HIV2> getDrugResistAlgorithm(String name) {
		return dl.getDrugResistAlgorithm(name);
	}

	@Override
	public DrugResistanceAlgorithm<HIV2> getDrugResistAlgorithm(String family, String version) {
		return dl.getDrugResistAlgorithm(family, version);
	}
	
	@Override
	public Gene<HIV2> extractMutationGene(String mutText) {
		return dl.extractMutationGene(mutText);
	}

	@Override
	public Mutation<HIV2> parseMutationString(Gene<HIV2> defaultGene, String mutText) {
		return dl.parseMutationString(defaultGene, mutText);
	}

	@Override
	public Mutation<HIV2> parseMutationString(String mutText) {
		return dl.parseMutationString(mutText);
	}
	
	@Override
	public MutationSet<HIV2> newMutationSet(String formattedMuts) {
		return dl.newMutationSet(formattedMuts);
	}

	@Override
	public MutationSet<HIV2> newMutationSet(Collection<String> formattedMuts) {
		return dl.newMutationSet(formattedMuts); 
	}

	@Override
	public MutationSet<HIV2> newMutationSet(Gene<HIV2> defaultGene, String formattedMuts) {
		return dl.newMutationSet(defaultGene, formattedMuts);
	}

	@Override
	public MutationSet<HIV2>	newMutationSet(Gene<HIV2> defaultGene, Collection<String> formattedMuts) {
		return dl.newMutationSet(defaultGene, formattedMuts);
	}

	@Override
	public Map<DrugClass<HIV2>, MutationSet<HIV2>> getDrugResistMutations() {
		return dl.getDrugResistMutations();
	}
	
	@Override
	public Map<DrugClass<HIV2>, MutationSet<HIV2>> getSurveilDrugResistMutations() {
		return dl.getSurveilDrugResistMutations();
	}

	@Override
	public Map<DrugClass<HIV2>, MutationSet<HIV2>> getRxSelectedMutations() {
		return dl.getRxSelectedMutations();
	}
	
	@Override
	public MutationSet<HIV2> getApobecMutations() {
		return dl.getApobecMutations();
	}

	@Override
	public MutationSet<HIV2> getApobecDRMs() {
		return dl.getApobecDRMs();
	}

	@Override
	public Collection<MutationType<HIV2>> getMutationTypes() {
		return dl.getMutationTypes();
	}
	
	@Override
	public MutationType<HIV2> getMutationType(String mutTypeText) {
		return dl.getMutationType(mutTypeText);
	}

	@Override
	public Collection<MutationTypePair<HIV2>> getMutationTypePairs() {
		return dl.getMutationTypePairs();
	}
	
	@Override
	public AminoAcidPercents<HIV2> getAminoAcidPercents(Strain<HIV2> strain, String treatment, String subtype) {
		return dl.getAminoAcidPercents(strain, treatment, subtype);
	}

	@Override
	public CodonPercents<HIV2> getCodonPercents(Strain<HIV2> strain, String treatment, String subtype) {
		return dl.getCodonPercents(strain, treatment, subtype);
	}

	@Override
	public List<MutationPrevalence<HIV2>> getMutationPrevalence(GenePosition<HIV2> genePos) {
		return dl.getMutationPrevalence(genePos);
	}
	
	@Override
	public ConditionalComments<HIV2> getConditionalComments() {
		return dl.getConditionalComments();
	}
	
	@Override
	public List<String> getMainSubtypes(Strain<HIV2> strain) {
		return dl.getMainSubtypes(strain);
	}
	
	@Override
	public Map<Gene<HIV2>, Map<String, Integer[]>> getNumPatientsForAAPercents(Strain<HIV2> strain) {
		return dl.getNumPatientsForAAPercents(strain);
	}

	@Override
	public Collection<Genotype<HIV2>> getGenotypes() {
		return dl.getGenotypes();
	}
	
	@Override
	public Genotype<HIV2> getGenotype(String name) {
		return dl.getGenotype(name);
	}

	@Override
	public Genotype<HIV2> getGenotypeUnknown() {
		return dl.getGenotypeUnknown();
	}

	@Override
	public List<GenotypeReference<HIV2>> getGenotypeReferences() {
		return dl.getGenotypeReferences();
	}
	
	@Override
	public Genotyper<HIV2> getGenotyper() {
		return dl.getGenotyper();
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// HIV instance is a singleton
		return false;
	}

	@Override
	public AminoAcidPercents<HIV2> getMainAminoAcidPercents(Strain<HIV2> strain) {
		return getAminoAcidPercents(strain, "all", "all");
	}

	@Override
	public CodonPercents<HIV2> getMainCodonPercents(Strain<HIV2> strain) {
		return getCodonPercents(strain, "all", "all");
	}

	@Override
	public DrugResistanceAlgorithm<HIV2> getDefaultDrugResistAlgorithm() {
		return getLatestDrugResistAlgorithm("HIVDB");
	}

	@Override
	public AlignmentConfig<HIV2> getAlignmentConfig() {
		return dl.getAlignmentConfig();
	}

	@Override
	public SequenceReadsAssembler<HIV2> getSequenceReadsAssembler(Strain<HIV2> strain) {
		return dl.getSequenceReadsAssemblers().get(strain);
	}
}