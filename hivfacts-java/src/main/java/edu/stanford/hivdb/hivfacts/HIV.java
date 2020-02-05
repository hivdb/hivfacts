package edu.stanford.hivdb.hivfacts;

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
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

public class HIV implements Virus<HIV> {

	private static final String VIRUS_NAME = "HIV";
	
	private static final String HIV1_STRAINS_RESPATH = "strains_hiv1.json";
	private static final String HIV1_GENES_RESPATH = "genes_hiv1.json";
	private static final String HIV1_DRUG_CLASSES_RESPATH = "drug-classes_hiv1.json";
	private static final String HIV_DRUGS_RESPATH = "drugs.json";
	private static final String HIV1_DRMS_RESPATH = "drms_hiv1.json";
	private static final String HIV1_SDRMS_RESPATH = "sdrms_hiv1.json";
	private static final String HIV1_TSMS_RESPATH = "tsms_hiv1.json";
	private static final String HIV1_APOBECS_RESPATH = "apobecs/apobecs.json";
	private static final String HIV1_APOBEC_DRMS_RESPATH = "apobecs/apobec_drms.json";
	private static final String HIV1_AAPCNTS_RESPATH = "aapcnt/rx-%s_subtype-%s.json";
	private static final String HIV1_CODONPCNTS_RESPATH = "codonpcnt/rx-%s_subtype-%s.json";
	private static final String HIV_MUTTYPES_RESPATH = "mutation-types.json";
	private static final String HIV1_MUTTYPE_PAIRS_RESPATH = "mutation-type-pairs_hiv1.json";
	private static final String HIV1_MAIN_SUBTYPES_RESPATH = "main-subtypes_hiv1.json";
	private static final String HIV1_GENOTYPE_REFERENCES_RESPATH = "genotypes/genotype-references_hiv1.json";
	private static final String HIV1_GENOTYPES_RESPATH = "genotypes/genotypes_hiv1.json";
	private static final String HIV_ALGORITHMS_INDEXPATH = "algorithms/versions.json";
	private static final String HIV_ALGORITHMS_RESPATH = "algorithms/%s_%s.xml";
	private static final String HIV1_CONDCOMMENTS_RESPATH = "conditional-comments_hiv1.json";

	private static final Pattern HIV_MUTATION_PATTERN = Pattern.compile(
		"^\\s*" +
		"((?i:PR|RT|IN))?[:_-]?" +
		"([AC-IK-NP-TV-Y])?" +
		"(\\d{1,3})" +
		"([AC-IK-NP-TV-Z.*]+(?:[#_]?[AC-IK-NP-TV-Z.*]+)?|[id_#~-]|[iI]ns(?:ertion)?|[dD]el(?:etion)?)" +
		"(?::([ACGTRYMWSKBDHVN-]{3})?)?" +
		"\\s*$");
	
	static {
		Virus.registerInstance(new HIV());
	}
		
	public static HIV getInstance() {
		return Virus.getInstance(HIV.class);
	}

	protected static String loadResource(String resPath) {
		try (
			InputStream stream = HIV.class
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
	
	private transient Map<String, Strain<HIV>> strains;
	private transient Map<String, Gene<HIV>> genes;
	private transient Map<String, DrugClass<HIV>> drugClasses;
	private transient Map<String, Drug<HIV>> drugs;
	private transient Map<DrugClass<HIV>, MutationSet<HIV>> drugResistMutations;
	private transient Map<DrugClass<HIV>, MutationSet<HIV>> surveilDrugResistMuts;
	private transient Map<DrugClass<HIV>, MutationSet<HIV>> rxSelectedMutations;
	private transient MutationSet<HIV> apobecMutations;
	private transient MutationSet<HIV> apobecDRMs;
	private transient Map<String, AminoAcidPercents<HIV>> aminoAcidPcnts = new HashMap<>();
	private transient Map<String, CodonPercents<HIV>> codonPcnts = new HashMap<>();
	private transient Map<String, MutationType<HIV>> mutationTypes;
	private transient List<MutationTypePair<HIV>> mutationTypePairs;
	private transient Map<Strain<HIV>, List<String>> mainSubtypes;
	private transient Map<GenePosition<HIV>, List<MutationPrevalence<HIV>>> mutPrevalenceMap = new HashMap<>();
	private transient Map<Strain<HIV>, Map<Gene<HIV>, Map<String, Integer[]>>> allAAPcntsNumPatients = new HashMap<>();
	private transient Map<String, Genotype<HIV>> allGenotypes;
	private transient List<GenotypeReference<HIV>> allGenotypeReferences;
	private transient Genotyper<HIV> genotyper;
	private transient List<DrugResistanceAlgorithm<HIV>> drugResistAlgs;
	private transient Map<String, DrugResistanceAlgorithm<HIV>> drugResistAlgLookup;
	private transient ConditionalComments<HIV> condComments;
	
	private HIV() {
		registerSequenceValidator(new HIVDefaultSequenceValidator());
		registerMutationsValidator(new HIVDefaultMutationsValidator());
		registerSequenceReadsValidator(new HIVDefaultSequenceReadsValidator());
	}

	private MutationSet<HIV> loadMutationSetFromRes(String resPath, Strain<HIV> strain) {
		String raw = loadResource(resPath);
		return MutationSet.loadJson(raw, geneText -> getGene(strain.toString() + geneText));
	}
	
	private Map<DrugClass<HIV>, MutationSet<HIV>> loadMutationSetByDrugClassFromRes(String resPath, Strain<HIV> strain) {
		Map<DrugClass<HIV>, MutationSet<HIV>> mutationsMap = new LinkedHashMap<>();
		String raw = loadResource(resPath);
		
		Map<String, List<Map<String, ?>>> muts = Json.loads(
			raw, new TypeToken<Map<String, List<Map<String, ?>>>>(){});
		for (String drugClassText : muts.keySet()) {
			DrugClass<HIV> drugClass = getDrugClass(drugClassText);
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
		String raw = loadResource(HIV1_CONDCOMMENTS_RESPATH);
		this.condComments = new ConditionalComments<>(raw, this);
	}
	
	private void initMainSubtypes() {
		String raw = loadResource(HIV1_MAIN_SUBTYPES_RESPATH);
		Map<String, List<String>> subtypes = Json.loads(raw, new TypeToken<Map<String, List<String>>>() {});
		Map<Strain<HIV>, List<String>> mainSubtypes = new LinkedHashMap<>();
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
		String raw = loadResource(HIV1_MUTTYPE_PAIRS_RESPATH);
		mutationTypePairs = MutationTypePair.loadJson(raw, this);
	}

	private void initStrains() {
		String raw = loadResource(HIV1_STRAINS_RESPATH);
		this.strains = Strain.loadJson(raw, this);
	}
	
	private void initGenes() {
		String raw = loadResource(HIV1_GENES_RESPATH);
		this.genes = Gene.loadJson(raw, this);
	}
	
	private void initDrugClasses() {
		String raw = loadResource(HIV1_DRUG_CLASSES_RESPATH);
		this.drugClasses = DrugClass.loadJson(raw, this);
	}
	
	private void initDrugs() {
		String raw = loadResource(HIV_DRUGS_RESPATH);
		this.drugs = Drug.loadJson(raw, this);
	}
	
	private void initDrugResistAlgs() {
		String raw = loadResource(HIV_ALGORITHMS_INDEXPATH);
		Map<String, List<List<String>>> algs = Json.loads(
			raw, new TypeToken<Map<String, List<List<String>>>>(){}.getType());
		List<DrugResistanceAlgorithm<HIV>> algList = new ArrayList<>();
		Map<String, DrugResistanceAlgorithm<HIV>> algMap = new LinkedHashMap<>(); 
		for (String family : algs.keySet()) {
			for (List<String> algData : algs.get(family)) {
				String version = algData.get(0);
				String publishDate = algData.get(1);
				String name = String.format("%s_%s", family, version);
				Strain<HIV> strain = getStrain(algData.get(2));
				String xmlText = loadResource(String.format(HIV_ALGORITHMS_RESPATH, family, version));
				DrugResistanceAlgorithm<HIV> alg = new DrugResistanceAlgorithm<>(
					name, family, version, publishDate, strain, xmlText);
				algList.add(alg);
				algMap.put(name, alg);
				algMap.put(alg.getEnumCompatName(), alg);
			}
		}
		this.drugResistAlgs = Collections.unmodifiableList(algList);
		this.drugResistAlgLookup = Collections.unmodifiableMap(algMap);
	}
	
	private void initGenotypes() {
		String raw = loadResource(HIV1_GENOTYPES_RESPATH);
		this.allGenotypes = Genotype.loadJson(raw, this);
	}
	
	private void initGenotypeReferences() {
		String raw = loadResource(HIV1_GENOTYPE_REFERENCES_RESPATH);
		this.allGenotypeReferences = GenotypeReference.loadJson(raw, this);
	}
	
	private void initDrugResistMutations() {
		drugResistMutations = loadMutationSetByDrugClassFromRes(HIV1_DRMS_RESPATH, getStrain("HIV1"));
	}
	
	private void initSurveilDrugResistMuts() {
		surveilDrugResistMuts = loadMutationSetByDrugClassFromRes(HIV1_SDRMS_RESPATH, getStrain("HIV1"));
	}
	
	private void initApobecMutations() {
		apobecMutations = loadMutationSetFromRes(HIV1_APOBECS_RESPATH, getStrain("HIV1"));
	}
	
	private void initApobecDRMs() {
		apobecDRMs = loadMutationSetFromRes(HIV1_APOBEC_DRMS_RESPATH, getStrain("HIV1"));
	}
	
	private void initRxSelectedMutations() {
		this.rxSelectedMutations = loadMutationSetByDrugClassFromRes(HIV1_TSMS_RESPATH, getStrain("HIV1"));
	}
	
	@Override
	public String getName() {
		return VIRUS_NAME;
	}
	
	@Override
	public Collection<Strain<HIV>> getStrains() {
		if (strains == null) {
			initStrains();
		}
		return strains.values();
	}
	
	@Override
	public Strain<HIV> getStrain(String name) {
		if (strains == null) {
			initStrains();
		}
		return strains.get(name);
	}

	@Override
	public Collection<Gene<HIV>> getGenes(Strain<HIV> strain) {
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
	public Gene<HIV> getGene(String name) {
		if (genes == null) {
			initGenes();
		}
		return genes.get(name);
	}
	
	@Override
	public Collection<DrugClass<HIV>> getDrugClasses() {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return Sets.newLinkedHashSet(drugClasses.values());
	}
	
	@Override
	public Map<String, DrugClass<HIV>> getDrugClassSynonymMap() {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return drugClasses;
	}
	
	@Override
	public DrugClass<HIV> getDrugClass(String name) {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return drugClasses.get(name);
	}
	
	@Override
	public Collection<Drug<HIV>> getDrugs() {
		if (drugs == null) {
			initDrugs();
		}
		return Sets.newTreeSet(drugs.values());
	}
	
	@Override
	public Map<String, Drug<HIV>> getDrugSynonymMap() {
		if (drugs == null) {
			initDrugs();
		}
		return drugs;
	}

	@Override
	public Collection<DrugResistanceAlgorithm<HIV>> getDrugResistAlgorithms() {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return drugResistAlgs;
	}

	@Override
	public Collection<DrugResistanceAlgorithm<HIV>> getDrugResistAlgorithms(Collection<String> algorithmNames) {
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
	public DrugResistanceAlgorithm<HIV> getDrugResistAlgorithm(String name) {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return drugResistAlgLookup.get(name);
	}

	@Override
	public DrugResistanceAlgorithm<HIV> getDrugResistAlgorithm(String family, String version) {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return drugResistAlgLookup.get(String.format("%s_%s", family, version));
	}
	
	/**
	 * Extracts gene from mutText string
	 * @param mutText
	 * @return a Gene enum object
	 */
	@Override
	public Gene<HIV> extractMutationGene(String mutText) {
		Gene<HIV> gene = null;
		Matcher m = HIV_MUTATION_PATTERN.matcher(mutText);
		if (m.matches()) {
			try {
				// TODO: currently we only support parsing HIV1 mutations.
				//       need to design a new format for HIV2 mutations.
				gene = getGene("HIV1" + m.group(1).toUpperCase());
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
	public Mutation<HIV> parseMutationString(Gene<HIV> defaultGene, String mutText) {
		Matcher m = HIV_MUTATION_PATTERN.matcher(mutText);
		ConsensusMutation<HIV> mut = null;
		if (m.matches()) {
			Gene<HIV> gene;
			try {
				// TODO: currently we only support parsing HIV1 mutations.
				//       need to design a new format for HIV2 mutations.
				gene = getGene("HIV1" + m.group(1).toUpperCase());
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
	public Mutation<HIV> parseMutationString(String mutText) {
		return parseMutationString(null, mutText);
	}
	
	@Override
	public MutationSet<HIV> newMutationSet(String formattedMuts) {
		return newMutationSet(null, formattedMuts);
	}

	@Override
	public MutationSet<HIV> newMutationSet(Collection<String> formattedMuts) {
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
	public MutationSet<HIV> newMutationSet(Gene<HIV> defaultGene, String formattedMuts) {
		if (formattedMuts == null) {
			return new MutationSet<>();
		}
		return newMutationSet(
			defaultGene,
			Arrays.asList(formattedMuts.split("[\\s,;+\\.]+"))
		);
	}

	@Override
	public MutationSet<HIV>	newMutationSet(Gene<HIV> defaultGene, Collection<String> formattedMuts) {
		return MutationSet.parseString(
			defaultGene, formattedMuts, (gene, mStr) -> parseMutationString(gene, mStr));
	}

	@Override
	public Map<DrugClass<HIV>, MutationSet<HIV>> getDrugResistMutations() {
		if (drugResistMutations == null) {
			initDrugResistMutations();
		}
		return drugResistMutations;
	}
	
	@Override
	public Map<DrugClass<HIV>, MutationSet<HIV>> getSurveilDrugResistMutations() {
		if (surveilDrugResistMuts == null) {
			initSurveilDrugResistMuts();
		}
		return surveilDrugResistMuts;
	}

	@Override
	public Map<DrugClass<HIV>, MutationSet<HIV>> getRxSelectedMutations() {
		if (rxSelectedMutations == null) {
			initRxSelectedMutations();
		}
		return rxSelectedMutations;
	}
	
	@Override
	public MutationSet<HIV> getApobecMutations() {
		if (apobecMutations == null) {
			initApobecMutations();
		}
		return apobecMutations;
	}

	@Override
	public MutationSet<HIV> getApobecDRMs() {
		if (apobecDRMs == null) {
			initApobecDRMs();
		}
		return apobecDRMs;
	}

	@Override
	public Collection<MutationType<HIV>> getMutationTypes() {
		if (mutationTypes == null) {
			initMutationTypes();
		}
		return mutationTypes.values();
	}
	
	@Override
	public MutationType<HIV> getMutationType(String mutTypeText) {
		if (mutationTypes == null) {
			initMutationTypes();
		}
		return mutationTypes.get(mutTypeText);
	}

	@Override
	public Collection<MutationTypePair<HIV>> getMutationTypePairs() {
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
	public AminoAcidPercents<HIV> getAminoAcidPercents(Strain<HIV> strain, String treatment, String subtype) {
		String resourceName = String.format(HIV1_AAPCNTS_RESPATH, treatment, subtype);
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
	public CodonPercents<HIV> getCodonPercents(Strain<HIV> strain, String treatment, String subtype) {
		String resourceName = String.format(HIV1_CODONPCNTS_RESPATH, treatment, subtype);
		if (!codonPcnts.containsKey(resourceName)) {
			codonPcnts.put(resourceName, new CodonPercents<>(resourceName, this, strain));
			// Example of emptyInstance:
			// codonPcnts.put(resourceName, CodonPercents.newEmptyInstance());
		}
		return codonPcnts.get(resourceName);
	}

	@Override
	public List<MutationPrevalence<HIV>> getMutationPrevalence(GenePosition<HIV> genePos) {
		if (!mutPrevalenceMap.containsKey(genePos)) {
			mutPrevalenceMap.put(genePos, Virus.super.getMutationPrevalence(genePos));
		}
		return new ArrayList<>(mutPrevalenceMap.get(genePos));					
	}
	
	@Override
	public ConditionalComments<HIV> getConditionalComments() {
		if (condComments == null) {
			initCondComments();
		}
		return condComments;
	}
	
	@Override
	public List<String> getMainSubtypes(Strain<HIV> strain) {
		if (mainSubtypes == null) {
			initMainSubtypes();
		}
		return mainSubtypes.get(strain);
	}
	
	@Override
	public Map<Gene<HIV>, Map<String, Integer[]>> getNumPatientsForAAPercents(Strain<HIV> strain) {
		if (!allAAPcntsNumPatients.containsKey(strain)) {
			allAAPcntsNumPatients.put(strain, Virus.super.getNumPatientsForAAPercents(strain));
		}
		return allAAPcntsNumPatients.get(strain);
	}

	@Override
	public Collection<Genotype<HIV>> getGenotypes() {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.values();
	}
	
	@Override
	public Genotype<HIV> getGenotype(String name) {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.get(name);
	}

	@Override
	public Genotype<HIV> getGenotypeUnknown() {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.get("U");
	}

	@Override
	public List<GenotypeReference<HIV>> getGenotypeReferences() {
		if (allGenotypeReferences == null) {
			initGenotypeReferences();
		}
		return allGenotypeReferences;
	}
	
	@Override
	public Genotyper<HIV> getGenotyper() {
		if (genotyper == null) {
			genotyper = new Genotyper<>(this);
		}
		return genotyper;
	}
	
	@Override
	public Strain<HIV> getMainStrain() {
		return getStrain("HIV1");
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) { return true; }
		// HIV instance is a singleton
		return false;
	}

}