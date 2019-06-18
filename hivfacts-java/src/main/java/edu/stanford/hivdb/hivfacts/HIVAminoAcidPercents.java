package edu.stanford.hivdb.hivfacts;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
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
public class HIVAminoAcidPercents {
	
	final static protected Gson gson = new Gson();
	final static protected Map<String, HIVAminoAcidPercents> singletons = new HashMap<>();
	
	final protected List<HIVAminoAcidPercent> aminoAcidPcnts;
	final private Map<GenePosition, Map<Character, HIVAminoAcidPercent>> aminoAcidPcntMap = new HashMap<>();

	/**
	 * Get an HIVAminoAcidPercents instance
	 * 
	 * @param treatment "naive" or "art"
	 * @param subtype "all", "A", "B", "C", "D", "F", "G", "CRF01_AE", "CRF02_AG", "other"
	 */
	public static HIVAminoAcidPercents getInstance(String treatment, String subtype) {
		String resourceName = String.format("aapcnt/rx-%s_subtype-%s.json", treatment, subtype);
		if (!singletons.containsKey(resourceName)) {
			singletons.put(resourceName, new HIVAminoAcidPercents(resourceName));
		}
		return singletons.get(resourceName);
	}
	

	/**
	 * HIVAminoAcidPercents initializer
	 * 
	 * @param resourceName
	 */
	protected HIVAminoAcidPercents(String resourceName) {
		
		try (
			InputStream stream = this
				.getClass().getClassLoader()
				.getResourceAsStream(resourceName);
		) {
			String raw = IOUtils.toString(stream, StandardCharsets.UTF_8);
			aminoAcidPcnts = gson.fromJson(raw, new TypeToken<List<HIVAminoAcidPercent>>(){}.getType());
		} catch (IOException|NullPointerException e) {
			throw new ExceptionInInitializerError(
				String.format("Invalid resource name (%s)", resourceName)
			);
		}

		for (HIVAminoAcidPercent aaPcnt : aminoAcidPcnts) {
			GenePosition gp = aaPcnt.getGenePosition();
			aminoAcidPcntMap.putIfAbsent(gp, new LinkedHashMap<>());
			aminoAcidPcntMap.get(gp).put(aaPcnt.aa, aaPcnt);
		}
	}

	public List<HIVAminoAcidPercent> get() {
		// make a copy in case of any modification
		return new ArrayList<>(aminoAcidPcnts);
	}
	
	public List<HIVAminoAcidPercent> get(String gene) {
		return (aminoAcidPcnts
				.stream().filter(aap -> aap.gene.equals(gene))
				.collect(Collectors.toList()));
	}

	public List<HIVAminoAcidPercent> get(Enum<?> geneEnum) {
		String gene = geneEnum.toString();
		return get(gene);
	}

	public List<HIVAminoAcidPercent> get(String gene, int pos) {
		return new ArrayList<>(aminoAcidPcntMap.get(new GenePosition(gene, pos)).values());
	}

	public List<HIVAminoAcidPercent> get(Enum<?> geneEnum, int pos) {
		String gene = geneEnum.toString();
		return get(gene, pos);
	}

	public HIVAminoAcidPercent get(String gene, int pos, char aa) {
		return aminoAcidPcntMap.get(new GenePosition(gene, pos)).get(aa);
	}

	public HIVAminoAcidPercent get(Enum<?> geneEnum, int pos, char aa) {
		String gene = geneEnum.toString();
		return get(gene, pos, aa);
	}

	/**
	 * Returns the highest amino acid prevalence associated with each of
	 * the AA in a mixture.
	 *
	 * @param gene
	 * @param pos
	 * @param cons consensus at the position
	 * @param mixture
	 *
	 * @return Double highest amino acid prevalence
	 */
	public Double getHighestAAPercentValue(
		String gene, int pos, /* char cons,*/ String mixture
	) {
		Double pcntVal = 0.0;
		GenePosition gpos = new GenePosition(gene, pos);

		for (char aa : mixture.toCharArray()) {
			/* if (aa == cons || aa == '*') {
				// ignore consensus and stop codon
				continue;
			} */
			double aaPcntVal = aminoAcidPcntMap.get(gpos).get(aa).percent;
			pcntVal = Math.max(pcntVal, aaPcntVal);
		}
		return pcntVal;
	}
	
	public Double getHighestAAPercentValue(
		Enum<?> geneEnum, int pos, /* char cons,*/ String mixture
	) {
		String gene = geneEnum.toString();
		return getHighestAAPercentValue(gene, pos, /*cons,*/ mixture);
	}
	
	/**
	 * Returns true if the given mutation contains any unusual AA
	 * 
	 * @param gene
	 * @param pos
	 * @param aas
	 * @return true if contains unusual AA
	 */
	public Boolean containsUnusualAA(String gene, int pos, String aas) {
		GenePosition gpos = new GenePosition(gene, pos);
		for (char aa : aas.toCharArray()) {
			HIVAminoAcidPercent aaPcnt = aminoAcidPcntMap.get(gpos).get(aa);
			if (aaPcnt.isUnusual) {
				return true;
			}
		}
		return false;
	}
	
	public Boolean containsUnusualAA(Enum<?> geneEnum, int pos, String aas) {
		String gene = geneEnum.toString();
		return containsUnusualAA(gene, pos, aas);
	}

}
