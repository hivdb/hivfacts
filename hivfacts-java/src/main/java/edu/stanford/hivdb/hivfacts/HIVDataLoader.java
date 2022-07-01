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
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.commons.io.IOUtils;

import com.google.common.collect.Sets;
import com.google.gson.reflect.TypeToken;

import edu.stanford.hivdb.comments.ConditionalComments;
import edu.stanford.hivdb.drugresistance.algorithm.DrugResistanceAlgorithm;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.genotypes.Genotype;
import edu.stanford.hivdb.genotypes.GenotypeReference;
import edu.stanford.hivdb.genotypes.Genotyper;
import edu.stanford.hivdb.mutations.AAMutation;
import edu.stanford.hivdb.mutations.AminoAcidPercents;
import edu.stanford.hivdb.mutations.CodonPercents;
import edu.stanford.hivdb.mutations.CodonMutation;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationPrevalence;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.mutations.MutationTypePair;
import edu.stanford.hivdb.seqreads.SequenceReadsAssembler;
import edu.stanford.hivdb.sequences.AlignmentConfig;
import edu.stanford.hivdb.sequences.SequenceAssembler;
import edu.stanford.hivdb.utilities.AAUtils;
import edu.stanford.hivdb.utilities.AssertUtils;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.viruses.Virus;

public class HIVDataLoader<T extends Virus<T>> {

	private static final Pattern HIV_MUTATION_PATTERN = Pattern.compile(
		"^\\s*" +
		"(__ASI__)?((?i:CA|PR|RT|IN))?[:_-]?" +
		"([AC-IK-NP-TV-Y])?" +
		"(\\d{1,3})" +
		"([AC-IK-NP-TV-Zid.*]+(?:[#_]?[AC-IK-NP-TV-Z.*]+)?|[id_#~-]|[iI]ns(?:ertion)?|[dD]el(?:etion)?)" +
		"(?::([ACGTRYMWSKBDHVN-]{3})?)?" +
		"\\s*$");
	
	private static final Pattern NON_ASI_AA_PATTERN = Pattern.compile(
		"^([AC-IK-NP-TV-Z.*]+(?:[#_]?[AC-IK-NP-TV-Z.*]+)?|[id_#~-]|[iI]ns(?:ertion)?|[dD]el(?:etion)?)$"
	);

	protected static String loadResource(String resPath) {
		try (
			InputStream stream = HIVDataLoader.class
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

	private final T virus;
	private final String VIRUS_NAME;
	private final String MAIN_STRAIN;
	private final String STRAINS_RESPATH;
	private final String GENES_RESPATH;
	private final String DRUG_CLASSES_RESPATH;
	private final String DRUGS_RESPATH;
	private final String DRMS_RESPATH;
	private final String SDRMS_RESPATH;
	private final String TSMS_RESPATH;
	private final String APOBECS_RESPATH;
	private final String APOBEC_DRMS_RESPATH;
	private final String AAPCNTS_RESPATH;
	private final String CODONPCNTS_RESPATH;
	private final String MUTTYPES_RESPATH;
	private final String MUTTYPE_PAIRS_RESPATH;
	private final String MAIN_SUBTYPES_RESPATH;
	private final String GENOTYPE_REFERENCES_RESPATH;
	private final String GENOTYPES_RESPATH;
	private final String ALGORITHMS_INDEXPATH;
	private final String ALGORITHMS_RESPATH;
	private final String CONDCOMMENTS_RESPATH;
	private final String ALIGNCONFIG_RESPATH;
	private final String ASSEMBLYCONFIG_RESPATH;

	private transient Map<String, Strain<T>> strains;
	private transient Map<String, Gene<T>> genes;
	private transient Map<String, DrugClass<T>> drugClasses;
	private transient Map<String, Drug<T>> drugs;
	private transient Map<DrugClass<T>, MutationSet<T>> drugResistMutations;
	private transient Map<DrugClass<T>, MutationSet<T>> surveilDrugResistMuts;
	private transient Map<DrugClass<T>, MutationSet<T>> rxSelectedMutations;
	private transient MutationSet<T> apobecMutations;
	private transient MutationSet<T> apobecDRMs;
	private transient Map<String, AminoAcidPercents<T>> aminoAcidPcnts = new HashMap<>();
	private transient Map<String, CodonPercents<T>> codonPcnts = new HashMap<>();
	private transient Map<String, MutationType<T>> mutationTypes;
	private transient List<MutationTypePair<T>> mutationTypePairs;
	private transient Map<Strain<T>, List<String>> mainSubtypes;
	private transient Map<GenePosition<T>, List<MutationPrevalence<T>>> mutPrevalenceMap = new HashMap<>();
	private transient Map<Strain<T>, Map<Gene<T>, Map<String, Integer[]>>> allAAPcntsNumPatients = new HashMap<>();
	private transient Map<String, Genotype<T>> allGenotypes;
	private transient List<GenotypeReference<T>> allGenotypeReferences;
	private transient Genotyper<T> genotyper;
	private transient List<DrugResistanceAlgorithm<T>> drugResistAlgs;
	private transient Map<String, DrugResistanceAlgorithm<T>> drugResistAlgLookup;
	private transient ConditionalComments<T> condComments;
	private transient AlignmentConfig<T> alignmentConfig;
	private transient Map<Strain<T>, SequenceReadsAssembler<T>> sequenceReadsAssemblers;
	private transient Map<Strain<T>, SequenceAssembler<T>> sequenceAssemblers;
	
	public HIVDataLoader(
		T virus,
		final String VIRUS_NAME,
		final String MAIN_STRAIN,
		final String STRAINS_RESPATH,
		final String GENES_RESPATH,
		final String DRUG_CLASSES_RESPATH,
		final String DRUGS_RESPATH,
		final String DRMS_RESPATH,
		final String SDRMS_RESPATH,
		final String TSMS_RESPATH,
		final String APOBECS_RESPATH,
		final String APOBEC_DRMS_RESPATH,
		final String AAPCNTS_RESPATH,
		final String CODONPCNTS_RESPATH,
		final String MUTTYPES_RESPATH,
		final String MUTTYPE_PAIRS_RESPATH,
		final String MAIN_SUBTYPES_RESPATH,
		final String GENOTYPE_REFERENCES_RESPATH,
		final String GENOTYPES_RESPATH,
		final String ALGORITHMS_INDEXPATH,
		final String ALGORITHMS_RESPATH,
		final String CONDCOMMENTS_RESPATH,
		final String ALIGNCONFIG_RESPATH,
		final String ASSEMBLYCONFIG_RESPATH
	) {
		this.virus = virus;
		this.VIRUS_NAME = VIRUS_NAME;
		this.MAIN_STRAIN = MAIN_STRAIN;
		this.STRAINS_RESPATH = STRAINS_RESPATH;
		this.GENES_RESPATH = GENES_RESPATH;
		this.DRUG_CLASSES_RESPATH = DRUG_CLASSES_RESPATH;
		this.DRUGS_RESPATH = DRUGS_RESPATH;
		this.DRMS_RESPATH = DRMS_RESPATH;
		this.SDRMS_RESPATH = SDRMS_RESPATH;
		this.TSMS_RESPATH = TSMS_RESPATH;
		this.APOBECS_RESPATH = APOBECS_RESPATH;
		this.APOBEC_DRMS_RESPATH = APOBEC_DRMS_RESPATH;
		this.AAPCNTS_RESPATH = AAPCNTS_RESPATH;
		this.CODONPCNTS_RESPATH = CODONPCNTS_RESPATH;
		this.MUTTYPES_RESPATH = MUTTYPES_RESPATH;
		this.MUTTYPE_PAIRS_RESPATH = MUTTYPE_PAIRS_RESPATH;
		this.MAIN_SUBTYPES_RESPATH = MAIN_SUBTYPES_RESPATH;
		this.GENOTYPE_REFERENCES_RESPATH = GENOTYPE_REFERENCES_RESPATH;
		this.GENOTYPES_RESPATH = GENOTYPES_RESPATH;
		this.ALGORITHMS_INDEXPATH = ALGORITHMS_INDEXPATH;
		this.ALGORITHMS_RESPATH = ALGORITHMS_RESPATH;
		this.CONDCOMMENTS_RESPATH = CONDCOMMENTS_RESPATH;
		this.ALIGNCONFIG_RESPATH = ALIGNCONFIG_RESPATH;
		this.ASSEMBLYCONFIG_RESPATH = ASSEMBLYCONFIG_RESPATH;
	}
	
	private MutationSet<T> loadMutationSetFromRes(String resPath, Collection<Strain<T>> strains) {
		String raw = loadResource(resPath);
		return (
			strains.stream()
			.map(strain -> MutationSet.loadJson(raw, geneText -> strain.getGene(geneText)))
			.reduce(new MutationSet<>(), (acc, other) -> acc.mergesWith(other))
		);
	}
	
	private Map<DrugClass<T>, MutationSet<T>> loadMutationSetByDrugClassFromRes(String resPath, Collection<Strain<T>> strains) {
		Map<DrugClass<T>, MutationSet<T>> mutationsMap = new LinkedHashMap<>();
		String raw = loadResource(resPath);
		
		Map<String, List<Map<String, ?>>> muts = Json.loads(
			raw, new TypeToken<Map<String, List<Map<String, ?>>>>(){});
		for (String drugClassText : muts.keySet()) {
			DrugClass<T> drugClass = getDrugClass(drugClassText);
			mutationsMap.put(
				drugClass,
				strains.stream()
				.map(strain -> MutationSet.loadJsonMap(
					muts.get(drugClassText),
					geneText -> strain.getGene(geneText)
				))
				.reduce(new MutationSet<>(), (acc, other) -> acc.mergesWith(other))
			);
		}
		return Collections.unmodifiableMap(mutationsMap);
	}
	
	private void initCondComments() {
		String raw = loadResource(CONDCOMMENTS_RESPATH);
		this.condComments = new ConditionalComments<>(raw, virus);
	}
	
	private void initMainSubtypes() {
		String raw = loadResource(MAIN_SUBTYPES_RESPATH);
		Map<String, List<String>> subtypes = Json.loads(raw, new TypeToken<Map<String, List<String>>>() {});
		Map<Strain<T>, List<String>> mainSubtypes = new LinkedHashMap<>();
		for (Map.Entry<String, List<String>> entry : subtypes.entrySet()) {
			mainSubtypes.put(
				getStrain(entry.getKey()),
				Collections.unmodifiableList(entry.getValue()));
		}
		this.mainSubtypes = Collections.unmodifiableMap(mainSubtypes);
	}
	
	private void initMutationTypes() {
		String raw = loadResource(MUTTYPES_RESPATH);
		mutationTypes = MutationType.loadJson(raw, virus);
	}

	private void initMutationTypePairs() {
		String raw = loadResource(MUTTYPE_PAIRS_RESPATH);
		mutationTypePairs = MutationTypePair.loadJson(raw, virus);
	}

	private void initStrains() {
		String raw = loadResource(STRAINS_RESPATH);
		this.strains = Strain.loadJson(raw, virus);
	}
	
	private void initGenes() {
		String raw = loadResource(GENES_RESPATH);
		this.genes = Gene.loadJson(raw, virus);
	}
	
	private void initDrugClasses() {
		String raw = loadResource(DRUG_CLASSES_RESPATH);
		this.drugClasses = DrugClass.loadJson(raw, virus);
	}
	
	private void initDrugs() {
		String raw = loadResource(DRUGS_RESPATH);
		this.drugs = Drug.loadJson(raw, virus);
	}
	
	private void initDrugResistAlgs() {
		String raw = loadResource(ALGORITHMS_INDEXPATH);
		Map<String, List<List<String>>> algs = Json.loads(
			raw, new TypeToken<Map<String, List<List<String>>>>(){}.getType());
		List<DrugResistanceAlgorithm<T>> algList = new ArrayList<>();
		Map<String, DrugResistanceAlgorithm<T>> algMap = new LinkedHashMap<>(); 
		for (String family : algs.keySet()) {
			for (List<String> algData : algs.get(family)) {
				String version = algData.get(0);
				String publishDate = algData.get(1);
				String name = String.format("%s_%s", family, version);
				String xmlText = loadResource(String.format(ALGORITHMS_RESPATH, family, version));
				DrugResistanceAlgorithm<T> alg = new DrugResistanceAlgorithm<>(
					name, family, version, publishDate, virus, xmlText);
				algList.add(alg);
				algMap.put(name, alg);
				algMap.put(alg.getEnumCompatName(), alg);
			}
		}
		this.drugResistAlgs = Collections.unmodifiableList(algList);
		this.drugResistAlgLookup = Collections.unmodifiableMap(algMap);
	}
	
	private void initGenotypes() {
		String raw = loadResource(GENOTYPES_RESPATH);
		this.allGenotypes = Genotype.loadJson(raw, virus);
	}
	
	private void initGenotypeReferences() {
		String raw = loadResource(GENOTYPE_REFERENCES_RESPATH);
		this.allGenotypeReferences = GenotypeReference.loadJson(raw, virus);
	}
	
	private void initDrugResistMutations() {
		drugResistMutations = loadMutationSetByDrugClassFromRes(DRMS_RESPATH, getStrains());
	}
	
	private void initSurveilDrugResistMuts() {
		surveilDrugResistMuts = loadMutationSetByDrugClassFromRes(SDRMS_RESPATH, getStrains());
	}
	
	private void initApobecMutations() {
		apobecMutations = loadMutationSetFromRes(APOBECS_RESPATH, getStrains());
	}
	
	private void initApobecDRMs() {
		apobecDRMs = loadMutationSetFromRes(APOBEC_DRMS_RESPATH, getStrains());
	}
	
	private void initRxSelectedMutations() {
		this.rxSelectedMutations = loadMutationSetByDrugClassFromRes(TSMS_RESPATH, getStrains());
	}
	
	public String getName() {
		return VIRUS_NAME;
	}
	
	public Strain<T> getMainStrain() {
		return getStrain(MAIN_STRAIN);
	}
	
	public Collection<Strain<T>> getStrains() {
		if (strains == null) {
			initStrains();
		}
		return strains.values();
	}
	
	
	public Strain<T> getStrain(String name) {
		if (strains == null) {
			initStrains();
		}
		return AssertUtils.notNull(
			strains.get(name),
			"Strain \"%s\" not found", name
		);
	}

	
	public Collection<Gene<T>> getGenes(Strain<T> strain) {
		if (genes == null) {
			initGenes();
		}
		return (
			genes.values()
			.stream()
			.distinct()
			.filter(gene -> gene.getStrain() == strain)
			.collect(Collectors.toList())
		);
	}
	
	
	public Gene<T> getGene(String name) {
		if (genes == null) {
			initGenes();
		}
		return AssertUtils.notNull(
			genes.get(name),
			"Gene \"%s\" not found", name
		);
	}
	
	
	public Collection<DrugClass<T>> getDrugClasses() {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return drugClasses.values().stream()
			.distinct()
			.collect(
				Collectors.toCollection(LinkedHashSet::new)
			);
	}
	
	
	public Map<String, DrugClass<T>> getDrugClassSynonymMap() {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return drugClasses;
	}
	
	
	public DrugClass<T> getDrugClass(String name) {
		if (drugClasses == null) {
			initDrugClasses();
		}
		return drugClasses.get(name);
	}
	
	
	public Collection<Drug<T>> getDrugs() {
		if (drugs == null) {
			initDrugs();
		}
		return Sets.newTreeSet(drugs.values());
	}
	
	
	public Map<String, Drug<T>> getDrugSynonymMap() {
		if (drugs == null) {
			initDrugs();
		}
		return drugs;
	}

	
	public Collection<DrugResistanceAlgorithm<T>> getDrugResistAlgorithms() {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return drugResistAlgs;
	}

	
	public Collection<DrugResistanceAlgorithm<T>> getDrugResistAlgorithms(Collection<String> algorithmNames) {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return (
			algorithmNames.stream()
			.map(name -> drugResistAlgLookup.get(name))
			.collect(Collectors.toList())
		);
	}
	
	
	
	public DrugResistanceAlgorithm<T> getDrugResistAlgorithm(String name) {
		if (drugResistAlgs == null) {
			initDrugResistAlgs();
		}
		return AssertUtils.notNull(
			drugResistAlgLookup.get(name),
			"Unable to locate algorithm %s", name
		);
	}

	
	public DrugResistanceAlgorithm<T> getDrugResistAlgorithm(String family, String version) {
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
	
	public Gene<T> extractMutationGene(String mutText) {
		Gene<T> gene = null;
		Matcher m = HIV_MUTATION_PATTERN.matcher(mutText);
		if (m.matches()) {
			try {
				gene = getGene(MAIN_STRAIN + m.group(2).toUpperCase());
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
	
	public Mutation<T> parseMutationString(Gene<T> defaultGene, String mutText) {
		Matcher m = HIV_MUTATION_PATTERN.matcher(mutText);
		Mutation<T> mut = null;
		if (m.matches()) {
			Gene<T> gene;
			boolean isASI = m.group(1) == null ? false : m.group(1).equals("__ASI__");
			try {
				gene = getGene(MAIN_STRAIN + m.group(2).toUpperCase());
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
			int pos = Integer.parseInt(m.group(4));
			String aas = m.group(5);
			if (!isASI && !NON_ASI_AA_PATTERN.matcher(aas).matches()) {
				throw new Mutation.InvalidMutationException(
					"Tried to parse mutation string using invalid parameters: " + mutText);
			}
			aas = AAUtils.normalizeAAs(aas);
			String triplet = m.group(6);
			if (triplet == null) triplet = "";
			if (isASI) {
				mut = new AAMutation<>(gene, pos, aas.toCharArray());
			}
			else {
				mut = new CodonMutation<>(gene, pos, aas, triplet, "", 0xff);
			}
		} else {
			throw new Mutation.InvalidMutationException(
				"Tried to parse mutation string using invalid parameters: " + mutText);
		}
		return mut;
	}

	
	public Mutation<T> parseMutationString(String mutText) {
		return parseMutationString(null, mutText);
	}
	
	
	public MutationSet<T> newMutationSet(String formattedMuts) {
		return newMutationSet(null, formattedMuts);
	}

	
	public MutationSet<T> newMutationSet(Collection<String> formattedMuts) {
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
	
	public MutationSet<T> newMutationSet(Gene<T> defaultGene, String formattedMuts) {
		if (formattedMuts == null) {
			return new MutationSet<>();
		}
		return newMutationSet(
			defaultGene,
			Arrays.asList(formattedMuts.split("[\\s,;+\\.]+"))
		);
	}

	
	public MutationSet<T>	newMutationSet(Gene<T> defaultGene, Collection<String> formattedMuts) {
		return MutationSet.parseString(
			defaultGene, formattedMuts, (gene, mStr) -> parseMutationString(gene, mStr));
	}

	
	public Map<DrugClass<T>, MutationSet<T>> getDrugResistMutations() {
		if (drugResistMutations == null) {
			initDrugResistMutations();
		}
		return drugResistMutations;
	}
	
	
	public Map<DrugClass<T>, MutationSet<T>> getSurveilDrugResistMutations() {
		if (surveilDrugResistMuts == null) {
			initSurveilDrugResistMuts();
		}
		return surveilDrugResistMuts;
	}

	
	public Map<DrugClass<T>, MutationSet<T>> getRxSelectedMutations() {
		if (rxSelectedMutations == null) {
			initRxSelectedMutations();
		}
		return rxSelectedMutations;
	}
	
	
	public MutationSet<T> getApobecMutations() {
		if (apobecMutations == null) {
			initApobecMutations();
		}
		return apobecMutations;
	}

	
	public MutationSet<T> getApobecDRMs() {
		if (apobecDRMs == null) {
			initApobecDRMs();
		}
		return apobecDRMs;
	}

	
	public Collection<MutationType<T>> getMutationTypes() {
		if (mutationTypes == null) {
			initMutationTypes();
		}
		return mutationTypes.values();
	}
	
	
	public MutationType<T> getMutationType(String mutTypeText) {
		if (mutationTypes == null) {
			initMutationTypes();
		}
		return mutationTypes.get(mutTypeText);
	}

	
	public Collection<MutationTypePair<T>> getMutationTypePairs() {
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
	
	public AminoAcidPercents<T> getAminoAcidPercents(Strain<T> strain, String treatment, String subtype) {
		String resourceName = String.format(AAPCNTS_RESPATH, treatment, subtype);
		String resourceKey = String.format("%s::%s", resourceName, strain.getName());
		if (!aminoAcidPcnts.containsKey(resourceKey)) {
			aminoAcidPcnts.put(resourceKey, new AminoAcidPercents<>(resourceName, virus, strain));
			// Example of empty Instance:
			// aminoAcidPcnts.put(resourceName, AminoAcidPercents.newEmptyInstance());
		}
		return aminoAcidPcnts.get(resourceKey);
	}

	/**
	 * Get a CodonPercents instance
	 *
	 * @param treatment "all", "naive" or "art"
	 * @param subtype "all", "A", "B", "C", "D", "F", "G", "CRF01_AE", "CRF02_AG"
	 */
	
	public CodonPercents<T> getCodonPercents(Strain<T> strain, String treatment, String subtype) {
		String resourceName = String.format(CODONPCNTS_RESPATH, treatment, subtype);
		if (!codonPcnts.containsKey(resourceName)) {
			codonPcnts.put(resourceName, new CodonPercents<>(resourceName, virus, strain));
			// Example of emptyInstance:
			// codonPcnts.put(resourceName, CodonPercents.newEmptyInstance());
		}
		return codonPcnts.get(resourceName);
	}

	
	public List<MutationPrevalence<T>> getMutationPrevalence(GenePosition<T> genePos) {
		if (!mutPrevalenceMap.containsKey(genePos)) {
			mutPrevalenceMap.put(genePos, virus.defaultGetMutationPrevalence(genePos));
		}
		return new ArrayList<>(mutPrevalenceMap.get(genePos));					
	}
	
	
	public ConditionalComments<T> getConditionalComments() {
		if (condComments == null) {
			initCondComments();
		}
		return condComments;
	}
	
	
	public List<String> getMainSubtypes(Strain<T> strain) {
		if (mainSubtypes == null) {
			initMainSubtypes();
		}
		return mainSubtypes.get(strain);
	}
	
	
	public Map<Gene<T>, Map<String, Integer[]>> getNumPatientsForAAPercents(Strain<T> strain) {
		if (!allAAPcntsNumPatients.containsKey(strain)) {
			allAAPcntsNumPatients.put(strain, virus.defaultGetNumPatientsForAAPercents(strain));
		}
		return allAAPcntsNumPatients.get(strain);
	}

	
	public Collection<Genotype<T>> getGenotypes() {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.values();
	}
	
	
	public Genotype<T> getGenotype(String name) {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.get(name);
	}

	
	public Genotype<T> getGenotypeUnknown() {
		if (allGenotypes == null) {
			initGenotypes();
		}
		return allGenotypes.get("U");
	}

	
	public List<GenotypeReference<T>> getGenotypeReferences() {
		if (allGenotypeReferences == null) {
			initGenotypeReferences();
		}
		return allGenotypeReferences;
	}
	
	
	public Genotyper<T> getGenotyper() {
		if (genotyper == null) {
			genotyper = new Genotyper<>(virus);
		}
		return genotyper;
	}

	public AlignmentConfig<T> getAlignmentConfig() {
		if (alignmentConfig == null) {
			String raw = loadResource(ALIGNCONFIG_RESPATH);
			alignmentConfig = AlignmentConfig.loadJson(raw, virus);
		}
		return alignmentConfig;
	}
	
	public Map<Strain<T>, SequenceReadsAssembler<T>> getSequenceReadsAssemblers() {
		if (sequenceReadsAssemblers == null) {
			String raw = loadResource(ASSEMBLYCONFIG_RESPATH);
			sequenceReadsAssemblers = SequenceReadsAssembler.loadJson(raw, virus);
		}
		return sequenceReadsAssemblers;
	}
	
	public Map<Strain<T>, SequenceAssembler<T>> getSequenceAssemblers() {
		if (sequenceAssemblers == null) {
			String raw = loadResource(ASSEMBLYCONFIG_RESPATH);
			Map<Strain<T>, SequenceAssembler<T>> sequenceAssemblers = SequenceAssembler.loadJson(
				raw,
				virus);
			this.sequenceAssemblers = sequenceAssemblers;
		}
		return sequenceAssemblers;
	}
	
}
