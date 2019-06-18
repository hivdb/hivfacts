package edu.stanford.hivdb.hivfacts;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.io.IOUtils;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;


/**
 * There are two public methods: getHighestMutPrevalence, unusualMutations
 *
 */
public class HIVCodonPercents {
	
	final static protected Gson gson = new Gson();
	final static protected Map<String, HIVCodonPercents> singletons = new HashMap<>();
	
	final protected List<HIVCodonPercent> codonPcnts;
	final private Map<GenePosition, Map<String, HIVCodonPercent>> codonPcntMap = new HashMap<>();

	/**
	 * Get an HIVCodonPercents instance
	 * 
	 * @param treatment "naive" or "art"
	 * @param subtype "all", "A", "B", "C", "D", "F", "G", "CRF01_AE", "CRF02_AG"
	 */
	public static HIVCodonPercents getInstance(String treatment, String subtype) {
		String resourceName = String.format("codonpcnt/rx-%s_subtype-%s.json", treatment, subtype);
		if (!singletons.containsKey(resourceName)) {
			singletons.put(resourceName, new HIVCodonPercents(resourceName));
		}
		return singletons.get(resourceName);
	}
	

	/**
	 * HIVAminoAcidPercents initializer
	 * 
	 * @param resourceName
	 */
	protected HIVCodonPercents(String resourceName) {
		
		try (
			InputStream stream = this
				.getClass().getClassLoader()
				.getResourceAsStream(resourceName);
		) {
			String raw = IOUtils.toString(stream, StandardCharsets.UTF_8);
			codonPcnts = gson.fromJson(raw, new TypeToken<List<HIVCodonPercent>>(){}.getType());
		} catch (IOException|NullPointerException e) {
			throw new ExceptionInInitializerError(
				String.format("Invalid resource name (%s)", resourceName)
			);
		}

		for (HIVCodonPercent cdPcnt : codonPcnts) {
			GenePosition gp = cdPcnt.getGenePosition();
			codonPcntMap.putIfAbsent(gp, new LinkedHashMap<>());
			codonPcntMap.get(gp).put(cdPcnt.codon, cdPcnt);
		}
	}

	public List<HIVCodonPercent> get() {
		// make a copy in case of any modification
		return new ArrayList<>(codonPcnts);
	}
	
	public List<HIVCodonPercent> get(String gene) {
		return (codonPcnts
				.stream().filter(cdp -> cdp.gene.equals(gene))
				.collect(Collectors.toList()));
	}

	public List<HIVCodonPercent> get(Enum<?> geneEnum) {
		String gene = geneEnum.toString();
		return get(gene);
	}

	public List<HIVCodonPercent> get(String gene, int pos) {
		return new ArrayList<>(codonPcntMap.getOrDefault(new GenePosition(gene, pos), Collections.emptyMap()).values());
	}

	public List<HIVCodonPercent> get(Enum<?> geneEnum, int pos) {
		String gene = geneEnum.toString();
		return get(gene, pos);
	}

	public HIVCodonPercent get(String gene, int pos, String codon) {
		Map<String, HIVCodonPercent> posCodons =
			codonPcntMap.getOrDefault(new GenePosition(gene, pos), Collections.emptyMap());
		if (posCodons.containsKey(codon)) {
			return posCodons.get(codon);
		}
		else if (posCodons.isEmpty()) {
			return null;
		}
		else {
			int total = posCodons.values().iterator().next().total;
			return new HIVCodonPercent(gene, pos, codon, 'X', .0, 0, total);
		}
	}

	public HIVCodonPercent get(Enum<?> geneEnum, int pos, String codon) {
		String gene = geneEnum.toString();
		return get(gene, pos, codon);
	}

	/**
	 * Returns the highest codon prevalence associated with each of
	 * the codon in a mixture.
	 *
	 * @param gene
	 * @param pos
	 * @param codonMixture
	 *
	 * @return Double highest amino acid prevalence
	 */
	public Double getHighestAAPercentValue(
		String gene, int pos, String[] codonMixture
	) {
		Double pcntVal = 0.0;

		for (String cd : codonMixture) {
			double cdPcntVal = get(gene, pos, cd).percent;
			pcntVal = Math.max(pcntVal, cdPcntVal);
		}
		return pcntVal;
	}
	
	public Double getHighestAAPercentValue(
		Enum<?> geneEnum, int pos, String[] codonMixture
	) {
		String gene = geneEnum.toString();
		return getHighestAAPercentValue(gene, pos, codonMixture);
	}
	
}
