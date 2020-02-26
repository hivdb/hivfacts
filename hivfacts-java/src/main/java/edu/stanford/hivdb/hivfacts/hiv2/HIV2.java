package edu.stanford.hivdb.hivfacts.hiv2;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.commons.io.IOUtils;

import com.google.common.collect.Sets;
import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.comments.ConditionalComments;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.drugs.DrugResistanceAlgorithm;
import edu.stanford.hivdb.genotypes.Genotype;
import edu.stanford.hivdb.genotypes.GenotypeReference;
import edu.stanford.hivdb.genotypes.Genotyper;
import edu.stanford.hivdb.mutations.AminoAcidPercents;
import edu.stanford.hivdb.mutations.CodonPercents;
import edu.stanford.hivdb.mutations.ConsensusMutation;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationPrevalence;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.mutations.MutationTypePair;
import edu.stanford.hivdb.utilities.AAUtils;
import edu.stanford.hivdb.utilities.AssertUtils;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

public class HIV2 implements Virus<HIV2> {

	private static final String VIRUS_NAME = "HIV2";
	
	private static final String HIV2_STRAINS_RESPATH = "strains_hiv2.json";
	private static final String HIV2_GENES_RESPATH = "genes_hiv2.json";
	private static final String HIV2_DRUG_CLASSES_RESPATH = "drug-classes_hiv2.json";
	private static final String HIV_DRUGS_RESPATH = "drugs.json";
	private static final String HIV2_DRMS_RESPATH = "drms_hiv2.json";
	private static final String HIV2_SDRMS_RESPATH = "sdrms_hiv2.json";
	private static final String HIV2_TSMS_RESPATH = "tsms_hiv2.json";
	private static final String HIV2_APOBECS_RESPATH = "apobecs/apobecs_hiv2.json";
	private static final String HIV2_APOBEC_DRMS_RESPATH = "apobecs/apobec_drms_hiv2.json";
	private static final String HIV2_AAPCNTS_RESPATH = "aapcnt/rx-%s_subtype-%s.json";
	private static final String HIV2_CODONPCNTS_RESPATH = "codonpcnt/rx-%s_subtype-%s.json";
	private static final String HIV_MUTTYPES_RESPATH = "mutation-types.json";
	private static final String HIV2_MUTTYPE_PAIRS_RESPATH = "mutation-type-pairs_hiv2.json";
	private static final String HIV2_MAIN_SUBTYPES_RESPATH = "main-subtypes_hiv2.json";
	private static final String HIV2_GENOTYPE_REFERENCES_RESPATH = "genotypes/genotype-references_hiv2.json";
	private static final String HIV2_GENOTYPES_RESPATH = "genotypes/genotypes_hiv2.json";
	private static final String HIV2_ALGORITHMS_INDEXPATH = "algorithms/versions_hiv2.json";
	private static final String HIV2_ALGORITHMS_RESPATH = "algorithms/%s-HIV2_%s.xml";
	private static final String HIV2_CONDCOMMENTS_RESPATH = "conditional-comments_hiv2.json";

	private static final Pattern HIV_MUTATION_PATTERN = Pattern.compile(
		"^\\s*" +
		"((?i:PR|RT|IN))?[:_-]?" +
		"([AC-IK-NP-TV-Y])?" +
		"(\\d{1,3})" +
		"([AC-IK-NP-TV-Z.*]+(?:[#_]?[AC-IK-NP-TV-Z.*]+)?|[id_#~-]|[iI]ns(?:ertion)?|[dD]el(?:etion)?)" +
		"(?::([ACGTRYMWSKBDHVN-]{3})?)?" +
		"\\s*$");
	
	static {
		Virus.registerInstance(new HIV2());
	}
		
	public static HIV2 getInstance() {
		return Virus.getInstance(HIV2.class);
	}

	protected static String loadResource(String resPath) {
		try (
			InputStream stream = HIV2.class
				.getClassLoader()
				.getResourceAsStream(resPath);
		) {
			return IOUtils.toString(stream, StandardCharsets.UTF_8);
		} catch (IOException|NullPointerException e) {
			throw new ExceptionInInitializerError(
				String.format("Invalid resource name (%s)", resPath)
			);
		}
	}
	
	private transient Map<String, Strain<HIV2>> strains;
	private transient Map<String, Gene<HIV2>> genes;
	private transient Map<String, DrugClass<HIV2>> drugClasses;
	private transient Map<String, Drug<HIV2>> drugs;
	private transient Map<DrugClass<HIV2>, MutationSet<HIV2>> drugResistMutations;
	private transient Map<DrugClass<HIV2>, MutationSet<HIV2>> surveilDrugResistMuts;
	private transient Map<DrugClass<HIV2>, MutationSet<HIV2>> rxSelectedMutations;
	private transient MutationSet<HIV2> apobecMutations;
	private transient MutationSet<HIV2> apobecDRMs;
	private transient Map<String, AminoAcidPercents<HIV2>> aminoAcidPcnts = new HashMap<>();
	private transient Map<String, CodonPercents<HIV2>> codonPcnts = new HashMap<>();
	private transient Map<String, MutationType<HIV2>> mutationTypes;
	private transient List<MutationTypePair<HIV2>> mutationTypePairs;
	private transient Map<Strain<HIV2>, List<String>> mainSubtypes;
	private transient Map<GenePosition<HIV2>, List<MutationPrevalence<HIV2>>> mutPrevalenceMap = new HashMap<>();
	private transient Map<Strain<HIV2>, Map<Gene<HIV2>, Map<String, Integer[]>>> allAAPcntsNumPatients = new HashMap<>();
	private transient Map<String, Genotype<HIV2>> allGenotypes;
	private transient List<GenotypeReference<HIV2>> allGenotypeReferences;
	private transient Genotyper<HIV2> genotyper;
	private transient List<DrugResistanceAlgorithm<HIV2>> drugResistAlgs;
	private transient Map<String, DrugResistanceAlgorithm<HIV2>> drugResistAlgLookup;
	private transient ConditionalComments<HIV2> condComments;
	
	private HIV2() {
		registerSequenceValidator(new HIV2DefaultSequenceValidator());
		registerMutationsValidator(new HIV2DefaultMutationsValidator());
		registerSequenceReadsValidator(new HIV2DefaultSequenceReadsValidator());
	}

	private MutationSet<HIV2> loadMutationSetFromRes(String resPath, Strain<HIV2> strain) {
		String raw = loadResource(resPath);
		return MutationSet.loadJson(raw, geneText -> getGene(strain.toString() + geneText));
	}
	
	private Map<DrugClass<HIV2>, MutationSet<HIV2>> loadMutationSetByDrugClassFromRes(String resPath, Strain<HIV2> strain) {
		Map<DrugClass<HIV2>, MutationSet<HIV2>> mutationsMap = new LinkedHashMap<>();
		String raw = loadResource(resPath);
		
		Map<String, List<Map<String, ?>>> muts = Json.loads(
			raw, new TypeToken<Map<String, List<Map<String, ?>>>>(){});
		for (String drugClassText : muts.keySet()) {
			DrugClass<HIV2> drugClass = getDrugClass(drugClassText);
			mutationsMap.put(
				drugClass,
				MutationSet.loadJsonMap(
					muts.get(drugClassText),
					geneText -> strain.getGene(geneText)
				)
			);
		}
		return Collections.unmodifiableMap(mutationsMap);
	}
	
	private void initCondComments() {
		String raw = loadResource(HIV2_CONDCOMMENTS_RESPATH);
		this.condComments = new ConditionalComments<>(raw, this);
	}
	
	private void initMainSubtypes() {
		String raw = loadResource(HIV2_MAIN_SUBTYPES_RESPATH);
		Map<String, List<String>> subtypes = Json.loads(raw, new TypeToken<Map<String, List<String>>>() {});
		Map<Strain<HIV2>, List<String>> mainSubtypes = new LinkedHashMap<>();
		for (Map.Entry<String, List<String>> entry : subtypes.entrySet()) {
			mainSubtypes.put(
				getStrain(entry.getKey()),
				Collections.unmodifiableList(entry.getValue()));
		}
		this.mainSubtypes = Collections.unmodifiableMap(mainSubtypes);
	}
	
	private void initMutationTypes() {
		String raw = loadResource(HIV_MUTTYPES_RESPATH);
		mutationTypes = MutationType.loadJson(raw, this);
	}

	private void initMutationTypePairs() {
		String raw = loadResource(HIV2_MUTTYPE_PAIRS_RESPATH);
		mutationTypePairs = MutationTypePair.loadJson(raw, this);
	}

	private void initStrains() {
		String raw = loadResource(HIV2_STRAINS_RESPATH);
		this.strains = Strain.loadJson(raw, this);
	}
	
	private void initGenes() {
		String raw = loadResource(HIV2_GENES_RESPATH);
		this.genes = Gene.loadJson(raw, this);
	}
	
	private void initDrugClasses() {
		String raw = loadResource(HIV2_DRUG_CLASSES_RESPATH);
		this.drugClasses = DrugClass.loadJson(raw, this);
	}
	
	private void initDrugs() {
		String raw = loadResource(HIV_DRUGS_RESPATH);
		this.drugs = Drug.loadJson(raw, this);
	}
	
	private void initDrugResistAlgs() {
		String raw = loadResource(HIV2_ALGORITHMS_INDEXPATH);
		Map<String, List<List<String>>> algs = Json.loads(
			raw, new TypeToken<Map<String, List<List<String>>>>(){}.getType());
		List<DrugResistanceAlgorithm<HIV2>> algList = new ArrayList<>();
		Map<String, DrugResistanceAlgorithm<HIV2>> algMap = new LinkedHashMap<>(); 
		for (String family : algs.keySet()) {
			for (List<String> algData : algs.get(family)) {
				String version = algData.get(0);
				String publishDate = algData.get(1);
				String name = String.format("%s_%s", family, version);
				String xmlText = loadResource(String.format(HIV2_ALGORITHMS_RESPATH, family, version));
				DrugResistanceAlgorithm<HIV2> alg = new DrugResistanceAlgorithm<>(
					name, family, version, publishDate, this, xmlText);
				algList.add(alg);
				algMap.put(name, alg);
				algMap.put(alg.getEnumCompatName(), alg);
			}
		}
		this.drugResistAlgs = Collections.unmodifiableList(algList);
		this.drugResistAlgLookup = Collections.unmodifiableMap(algMap);
	}
	
	private void initGenotypes() {
		String raw = loadResource(HIV2_GENOTYPES_RESPATH);
		this.allGenotypes = Genotype.loadJson(raw, this);
	}
	
	private void initGenotypeReferences() {
		String raw = loadResource(HIV2_GENOTYPE_REFERENCES_RESPATH);
		this.allGenotypeReferences = GenotypeReference.loadJson(raw, this);
	}
	
	private void initDrugResistMutations() {
		drugResistMutations = loadMutationSetByDrugClassFromRes(HIV2_DRMS_RESPATH, getStrain("HIV2A"));
	}
	
	private void initSurveilDrugResistMuts() {
		surveilDrugResistMuts = loadMutationSetByDrugClassFromRes(HIV2_SDRMS_RESPATH, getStrain("HIV2A"));
	}
	
	private void initApobecMutations() {
		apobecMutations = loadMutationSetFromRes(HIV2_APOBECS_RESPATH, getStrain("HIV2A"));
	}
	
	private void initApobecDRMs() {
		apobecDRMs = loadMutationSetFromRes(HIV2_APOBEC_DRMS_RESPATH, getStrain("HIV2A"));
	}
	
	private void initRxSelectedMutations() {
		this.rxSelectedMutations = loadMutationSetByDrugClassFromRes(HIV2_TSMS_RESPATH, getStrain("HIV2A"));
	}
	
	@Override
	public String getName() {
		return VIRUS_NAME;
	}
	
	@Override
	public Collection<Strain<HIV2>> getStrains() {
		if (strains == null) {
			initStrains();
		}
		return strains.values();
	}
	
	@Override
	public Strain<HIV2> getStrain(String name) {
		if (strains == null) {
			initStrains();
		}
		return AssertUtils.notNull(
			strains.get(name),
			"Strain \"%s\" not found", name
		);
	}

	@Override
	public Collection<Gene<HIV2>> getGenes(Strain<HIV2> strain) {
		if (genes == null) {
			initGenes();
		}
		return (
			genes.values()
			.stream()
			.filter(gene -> gene.getStrain() == strain)
			.collect(Collectors.toList())
		);
	}
	
	@Override
	public Gene<HIV2> getGene(String name) {
		if (genes == null) {
			initGenes();
		}
		return genes.get(name);
	}
	
	@Override
	public Collection<DrugClass<HIV2>> getDrugClasses() {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return Sets.newLinkedHashSet(drugClasses.values());
	}
	
	@Override
	public Map<String, DrugClass<HIV2>> getDrugClassSynonymMap() {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return drugClasses;
	}
	
	@Override
	public DrugClass<HIV2> getDrugClass(String name) {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return drugClasses.get(name);
	}
	
	@Override
	public Collection<Drug<HIV2>> getDrugs() {
		if (drugs == null) {
			initDrugs();
		}
		return Sets.newTreeSet(drugs.values());
	}
	
	@Override
	public Map<String, Drug<HIV2>> getDrugSynonymMap() {
		if (drugs == null) {
			initDrugs();
		}
		return drugs;
	}

	@Override
	public Collection<DrugResistanceAlgorithm<HIV2>> getDrugResistAlgorithms() {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return drugResistAlgs;
	}

	@Override
	public Collection<DrugResistanceAlgorithm<HIV2>> getDrugResistAlgorithms(Collection<String> algorithmNames) {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return (
			algorithmNames.stream()
			.map(name -> drugResistAlgLookup.get(name))
			.collect(Collectors.toList())
		);
	}
	
	
	@Override
	public DrugResistanceAlgorithm<HIV2> getDrugResistAlgorithm(String name) {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return AssertUtils.notNull(
			drugResistAlgLookup.get(name),
			"Unable to locate algorithm %s", name
		);
	}

	@Override
	public DrugResistanceAlgorithm<HIV2> getDrugResistAlgorithm(String family, String version) {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return AssertUtils.notNull(
			drugResistAlgLookup.get(String.format("%s_%s", family, version)),
			"Unable to locate algorithm %s_%s", family, version
		);
	}
	
	/**
	 * Extracts gene from mutText string
	 * @param mutText
	 * @return a Gene enum object
	 */
	@Override
	public Gene<HIV2> extractMutationGene(String mutText) {
		Gene<HIV2> gene = null;
		Matcher m = HIV_MUTATION_PATTERN.matcher(mutText);
		if (m.matches()) {
			try {
				// TODO: currently we only support parsing HIV2 mutations.
				//       need to design a new format for HIV2 mutations.
				gene = getGene("HIV2A" + m.group(1).toUpperCase());
			} catch (NullPointerException e) {
				throw new Mutation.InvalidMutationException(
					"Gene is not specified and also not found in the " +
					"given text: " + mutText + ". The correct format " +
					"for an input mutation string is, for example, " +
					"RT:215Y.", e);
			}
		}
		return gene;
	}

	/**
	 * Converts gene and mutText string into a Mutation object
	 * mutText may or may not have a preceding consensus
	 * @param gene, mutText
	 * @return a Mutation object
	 */
	@Override
	public Mutation<HIV2> parseMutationString(Gene<HIV2> defaultGene, String mutText) {
		Matcher m = HIV_MUTATION_PATTERN.matcher(mutText);
		ConsensusMutation<HIV2> mut = null;
		if (m.matches()) {
			Gene<HIV2> gene;
			try {
				// TODO: currently we only support parsing HIV2 mutations.
				//       need to design a new format for HIV2 mutations.
				gene = getGene("HIV2A" + m.group(1).toUpperCase());
			} catch (NullPointerException e) {
				if (defaultGene == null) {
					throw new Mutation.InvalidMutationException(
						"Gene is not specified and also not found in the " +
						"given text: " + mutText + ". The correct format " +
						"for an input mutation string is, for example, " +
						"RT:215Y.", e);
				}
				else {
					gene = defaultGene;
				}
			}
			int pos = Integer.parseInt(m.group(3));
			String aas = AAUtils.normalizeAAs(m.group(4));
			String triplet = m.group(5);
			if (triplet == null) triplet = "";
			mut = new ConsensusMutation<>(gene, pos, aas, triplet, "", 0xff);
		} else {
			throw new Mutation.InvalidMutationException(
				"Tried to parse mutation string using invalid parameters: " + mutText);
		}
		return mut;
	}

	@Override
	public Mutation<HIV2> parseMutationString(String mutText) {
		return parseMutationString(null, mutText);
	}
	
	@Override
	public MutationSet<HIV2> newMutationSet(String formattedMuts) {
		return newMutationSet(null, formattedMuts);
	}

	@Override
	public MutationSet<HIV2> newMutationSet(Collection<String> formattedMuts) {
		return newMutationSet(null, formattedMuts);
	}

	/**
	 * Parse the given string then create mutation set.
	 *
	 * Supported delimiters:
	 * 	- space ( )
	 * 	- tabulation (\t)
	 * 	- new line (\n)
	 * 	- carriage return (\r)
	 * 	- dot (.)
	 * 	- comma (,)
	 * 	- semicolon (;)
	 * 	- plus (+)
	 *
	 * Supported single mutation formats:
	 * 	- with consensus (P1X)
	 * 	- without consensus (52R)
	 * 	- lowercase indels (69i or 44d)
	 * 	- dash/underscore indels (69_XX or 44-)
	 * 	- "hivdb" indels (69#XX or 44~)
	 * 	- word indels (69Insertion, 44Deletion)
	 * 	- stop codon (122*)
	 *
	 * All duplicated mutations are removed before returning.
	 *
	 * @param defaultGene
	 * @param fomattedMuts
	 * @return A list of Mutation objects
	 */
	@Override
	public MutationSet<HIV2> newMutationSet(Gene<HIV2> defaultGene, String formattedMuts) {
		if (formattedMuts == null) {
			return new MutationSet<>();
		}
		return newMutationSet(
			defaultGene,
			Arrays.asList(formattedMuts.split("[\\s,;+\\.]+"))
		);
	}

	@Override
	public MutationSet<HIV2>	newMutationSet(Gene<HIV2> defaultGene, Collection<String> formattedMuts) {
		return MutationSet.parseString(
			defaultGene, formattedMuts, (gene, mStr) -> parseMutationString(gene, mStr));
	}

	@Override
	public Map<DrugClass<HIV2>, MutationSet<HIV2>> getDrugResistMutations() {
		if (drugResistMutations == null) {
			initDrugResistMutations();
		}
		return drugResistMutations;
	}
	
	@Override
	public Map<DrugClass<HIV2>, MutationSet<HIV2>> getSurveilDrugResistMutations() {
		if (surveilDrugResistMuts == null) {
			initSurveilDrugResistMuts();
		}
		return surveilDrugResistMuts;
	}

	@Override
	public Map<DrugClass<HIV2>, MutationSet<HIV2>> getRxSelectedMutations() {
		if (rxSelectedMutations == null) {
			initRxSelectedMutations();
		}
		return rxSelectedMutations;
	}
	
	@Override
	public MutationSet<HIV2> getApobecMutations() {
		if (apobecMutations == null) {
			initApobecMutations();
		}
		return apobecMutations;
	}

	@Override
	public MutationSet<HIV2> getApobecDRMs() {
		if (apobecDRMs == null) {
			initApobecDRMs();
		}
		return apobecDRMs;
	}

	@Override
	public Collection<MutationType<HIV2>> getMutationTypes() {
		if (mutationTypes == null) {
			initMutationTypes();
		}
		return mutationTypes.values();
	}
	
	@Override
	public MutationType<HIV2> getMutationType(String mutTypeText) {
		if (mutationTypes == null) {
			initMutationTypes();
		}
		return mutationTypes.get(mutTypeText);
	}

	@Override
	public Collection<MutationTypePair<HIV2>> getMutationTypePairs() {
		if (mutationTypePairs == null) {
			initMutationTypePairs();
		}
		return mutationTypePairs;
	}
	
	/**
	 * Get an AminoAcidPercents instance
	 *
	 * @param treatment "all", "naive" or "art"
	 * @param subtype "all", "A", "B", "C", "D", "F", "G", "CRF01_AE", "CRF02_AG", "other"
	 */
	@Override
	public AminoAcidPercents<HIV2> getAminoAcidPercents(Strain<HIV2> strain, String treatment, String subtype) {
		String resourceName = String.format(HIV2_AAPCNTS_RESPATH, treatment, subtype);
		if (!aminoAcidPcnts.containsKey(resourceName)) {
			aminoAcidPcnts.put(resourceName, new AminoAcidPercents<>(resourceName, this, strain));
			// Example of empty Instance:
			// aminoAcidPcnts.put(resourceName, AminoAcidPercents.newEmptyInstance());
		}
		return aminoAcidPcnts.get(resourceName);
	}

	/**
	 * Get a CodonPercents instance
	 *
	 * @param treatment "all", "naive" or "art"
	 * @param subtype "all", "A", "B", "C", "D", "F", "G", "CRF01_AE", "CRF02_AG"
	 */
	@Override
	public CodonPercents<HIV2> getCodonPercents(Strain<HIV2> strain, String treatment, String subtype) {
		String resourceName = String.format(HIV2_CODONPCNTS_RESPATH, treatment, subtype);
		if (!codonPcnts.containsKey(resourceName)) {
			codonPcnts.put(resourceName, new CodonPercents<>(resourceName, this, strain));
			// Example of emptyInstance:
			// codonPcnts.put(resourceName, CodonPercents.newEmptyInstance());
		}
		return codonPcnts.get(resourceName);
	}

	@Override
	public List<MutationPrevalence<HIV2>> getMutationPrevalence(GenePosition<HIV2> genePos) {
		if (!mutPrevalenceMap.containsKey(genePos)) {
			mutPrevalenceMap.put(genePos, Virus.super.getMutationPrevalence(genePos));
		}
		return new ArrayList<>(mutPrevalenceMap.get(genePos));					
	}
	
	@Override
	public ConditionalComments<HIV2> getConditionalComments() {
		if (condComments == null) {
			initCondComments();
		}
		return condComments;
	}
	
	@Override
	public List<String> getMainSubtypes(Strain<HIV2> strain) {
		if (mainSubtypes == null) {
			initMainSubtypes();
		}
		return mainSubtypes.get(strain);
	}
	
	@Override
	public Map<Gene<HIV2>, Map<String, Integer[]>> getNumPatientsForAAPercents(Strain<HIV2> strain) {
		if (!allAAPcntsNumPatients.containsKey(strain)) {
			allAAPcntsNumPatients.put(strain, Virus.super.getNumPatientsForAAPercents(strain));
		}
		return allAAPcntsNumPatients.get(strain);
	}

	@Override
	public Collection<Genotype<HIV2>> getGenotypes() {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.values();
	}
	
	@Override
	public Genotype<HIV2> getGenotype(String name) {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.get(name);
	}

	@Override
	public Genotype<HIV2> getGenotypeUnknown() {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.get("U");
	}

	@Override
	public List<GenotypeReference<HIV2>> getGenotypeReferences() {
		if (allGenotypeReferences == null) {
			initGenotypeReferences();
		}
		return allGenotypeReferences;
	}
	
	@Override
	public Genotyper<HIV2> getGenotyper() {
		if (genotyper == null) {
			genotyper = new Genotyper<>(this);
		}
		return genotyper;
	}
	
	@Override
	public Strain<HIV2> getMainStrain() {
		return getStrain("HIV2A");
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// HIV instance is a singleton
		return false;
	}

}