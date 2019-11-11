package edu.stanford.hivdb.hivfacts;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.IOUtils;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;

public class HIVAPOBECMutations {

	final static protected Gson gson = new Gson();
	final static protected HIVAPOBECMutations singleton;

	static {
		singleton = new HIVAPOBECMutations(
			"apobecs/apobecs.json",
			"apobecs/apobec_drms.json"
		);
	}

	final protected List<HIVAPOBECMutation> apobecs;
	final protected List<HIVAPOBECMutation> apobecDRMs;
	final private Map<GenePosition, Set<Character>> apobecMap;
	final private Map<GenePosition, Set<Character>> apobecDRMMap;

	public static HIVAPOBECMutations getInstance() {
		return singleton;
	}

	protected static List<HIVAPOBECMutation> loadAPOBECRes(String resPath) {
		try (
			InputStream stream = HIVAPOBECMutations.class
				.getClassLoader()
				.getResourceAsStream(resPath);
		) {
			String raw = IOUtils.toString(stream, StandardCharsets.UTF_8);
			return gson.fromJson(raw, new TypeToken<List<HIVAPOBECMutation>>(){}.getType());
		} catch (IOException|NullPointerException e) {
			throw new ExceptionInInitializerError(
				String.format("Invalid resource name (%s)", resPath)
			);
		}
	}

	protected static Map<GenePosition, Set<Character>>
			buildLookupMap(List<HIVAPOBECMutation> mutations) {
		Map<GenePosition, Set<Character>> lookup = new HashMap<>();
		for (HIVAPOBECMutation mut : mutations) {
			GenePosition gp = new GenePosition(mut.getGene(), mut.position);
			lookup.putIfAbsent(gp, new HashSet<>());
			lookup.get(gp).add(mut.aa);
		}
		return lookup;
	}

	/**
	 * HIVAminoAcidPercents initializer
	 *
	 * @param resourceName
	 */
	protected HIVAPOBECMutations(String resApobecsPath, String resApobecDRMsPath) {
		apobecs = loadAPOBECRes(resApobecsPath);
		apobecDRMs = loadAPOBECRes(resApobecDRMsPath);
		apobecMap = buildLookupMap(apobecs);
		apobecDRMMap = buildLookupMap(apobecDRMs);
	}

	public boolean isApobecMutation(Gene gene, int pos, char aa) {
		Set<Character> aas = apobecMap.get(new GenePosition(gene, pos));
		if (aas != null) {
			return aas.contains(aa);
		}
		return false;
	}

	public boolean isApobecDRM(Gene gene, int pos, char aa) {
		Set<Character> aas = apobecDRMMap.get(new GenePosition(gene, pos));
		if (aas != null) {
			return aas.contains(aa);
		}
		return false;
	}

	public List<HIVAPOBECMutation> getApobecMutations() { return apobecs; }
	public List<HIVAPOBECMutation> getApobecDRMs() { return apobecDRMs; }

}
