package edu.stanford.hivdb.aapcnt;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class HIVAminoAcidPercentsTest {
	
	private static enum Gene { PR, RT, IN }
	
	@Rule
	public ExpectedException expectedEx = ExpectedException.none();

	@Test
	public void testGetInstanceSuccess() {
		HIVAminoAcidPercents allAll01 = HIVAminoAcidPercents.getInstance("all", "All");
		HIVAminoAcidPercents allAll02 = HIVAminoAcidPercents.getInstance("all", "All");
		assertEquals("same singleton instance", allAll01, allAll02);
	}
	
	@Test
	public void testGetInstanceFailCase1() {
		expectedEx.expect(ExceptionInInitializerError.class);
		expectedEx.expectMessage("Invalid resource name (allE.json)");
		HIVAminoAcidPercents.getInstance("all", "E");
	}

	@Test
	public void testGetInstanceFailCase2() {
		expectedEx.expect(ExceptionInInitializerError.class);
		expectedEx.expectMessage("Invalid resource name (aaaaaaaAll.json)");
		HIVAminoAcidPercents.getInstance("aaaaaaa", "All");
	}
	
	@Test
	public void testGsonLoad() {
		HIVAminoAcidPercents allAll = HIVAminoAcidPercents.getInstance("all", "All");
		assertEquals((99 + 560 + 288) * 23, allAll.aminoAcidPcnts.size());
		HIVAminoAcidPercent mutPR1A = allAll.aminoAcidPcnts.get(0);
		assertEquals("PR", mutPR1A.gene);
		assertEquals(1, (int) mutPR1A.position);
		assertEquals('A', (char) mutPR1A.aa);
	}
	
	@Test
	public void testGetAll() {
		HIVAminoAcidPercents allAll = HIVAminoAcidPercents.getInstance("all", "All");
		List<HIVAminoAcidPercent> allAAPcnts = allAll.get();
		assertEquals(allAll.aminoAcidPcnts, allAAPcnts); // equal values
		assertFalse(allAll.aminoAcidPcnts == allAAPcnts); // different references
		assertEquals((99 + 560 + 288) * 23, allAAPcnts.size());
	}
	
	@Test
	public void testGetByGene() {
		HIVAminoAcidPercents allAll = HIVAminoAcidPercents.getInstance("all", "All");
		List<HIVAminoAcidPercent> gRTAAPcnts = allAll.get("RT");
		List<HIVAminoAcidPercent> gRTAAPcntsByEnum = allAll.get(Gene.RT);
		assertEquals(gRTAAPcnts, gRTAAPcntsByEnum); // equal values
		assertFalse(gRTAAPcnts == gRTAAPcntsByEnum); // different references
		assertEquals(560 * 23, gRTAAPcnts.size());
		int i = 0;
		for (int p = 1; p <= 560; p ++) {
			for (char aa : "ACDEFGHIKLMNPQRSTVWY_-*".toCharArray()) {
				HIVAminoAcidPercent mut = gRTAAPcnts.get(i ++);
				assertEquals("RT", mut.gene);
				assertEquals(p, (int) mut.position);
				assertEquals(aa, (char) mut.aa);
			}
		}
	}
	
	@Test
	public void testGetByGenePos() {
		HIVAminoAcidPercents allAll = HIVAminoAcidPercents.getInstance("all", "All");
		List<HIVAminoAcidPercent> gpIN263AAPcnts = allAll.get("IN", 263);
		List<HIVAminoAcidPercent> gpIN263AAPcntsByEnum = allAll.get(Gene.IN, 263);
		assertEquals(gpIN263AAPcnts, gpIN263AAPcntsByEnum); // equal value
		assertFalse(gpIN263AAPcnts == gpIN263AAPcntsByEnum); // different references
		assertEquals(23, gpIN263AAPcntsByEnum.size());
		int i = 0;
		for (char aa : "ACDEFGHIKLMNPQRSTVWY_-*".toCharArray()) {
			HIVAminoAcidPercent mut = gpIN263AAPcntsByEnum.get(i ++);
			assertEquals("IN", mut.gene);
			assertEquals(263, (int) mut.position);
			assertEquals(aa, (char) mut.aa);
		}
	}
	
	@Test
	public void testGetByMut() {
		HIVAminoAcidPercents allAll = HIVAminoAcidPercents.getInstance("all", "All");
		HIVAminoAcidPercent mutIN263R = allAll.get("IN", 263, 'R');
		HIVAminoAcidPercent mutIN263RByEnum = allAll.get(Gene.IN, 263, 'R');
		assertEquals(mutIN263R, mutIN263RByEnum); // equal value
		assertTrue(mutIN263R == mutIN263RByEnum); // same references
		assertEquals("IN", mutIN263R.gene);
		assertEquals(263, (int) mutIN263R.position);
		assertEquals('R', (char) mutIN263R.aa);
	}
	
	@Test
	public void testGetHighestPercentValue() {
		HIVAminoAcidPercents allAll = HIVAminoAcidPercents.getInstance("all", "All");
		double highestVal = allAll.getHighestAAPercentValue("PR", 23, "AHI");
		double highestValByEnum = allAll.getHighestAAPercentValue(Gene.PR, 23, "AHI");
		assertEquals(highestVal, highestValByEnum, 1e-18);
		double expectedHighestVal = .0;
		for (char aa : "AHI".toCharArray()) {
			double pcntVal = allAll.get("PR", 23, aa).percent;
			expectedHighestVal = Math.max(expectedHighestVal, pcntVal);
		}
		assertEquals(expectedHighestVal, highestVal, 1e-18);

		// These are intended to fail for every version update.
		// You must manually check and correct these numbers.
		assertEquals(0.00188219, highestVal, 1e-8);
		assertEquals(0.08008744, allAll.getHighestAAPercentValue("RT", 67, "N"), 1e-8);
		assertEquals(0.00825590, allAll.getHighestAAPercentValue("RT", 69, "KS"), 1e-8);
		assertEquals(0.04722477, allAll.getHighestAAPercentValue("PR", 82, "IA"), 1e-8);
		assertEquals(0.0, allAll.getHighestAAPercentValue(Gene.RT, 67, "W"), 1e-8);
		assertEquals(0.08008744, allAll.getHighestAAPercentValue("RT", 67, "N*"), 1e-8);
		assertEquals(0.0, allAll.getHighestAAPercentValue("RT", 67, "*"), 1e-8);
		assertEquals(0.0, allAll.getHighestAAPercentValue("IN", 1, ""), 1e-8);
	}
	
	@Test
	public void testContainsUnusualAA() {
		HIVAminoAcidPercents allAll = HIVAminoAcidPercents.getInstance("all", "All");
		assertTrue(allAll.containsUnusualAA("RT", 5, "I*"));
		assertTrue(allAll.containsUnusualAA(Gene.PR, 6, "W*")); // APOBEC stop codon
		assertFalse(allAll.containsUnusualAA("RT", 67, "N"));
		assertFalse(allAll.containsUnusualAA(Gene.RT, 69, "KS"));
		assertTrue(allAll.containsUnusualAA("PR", 82, "VIAD"));
		assertTrue(allAll.containsUnusualAA(Gene.RT, 67, "W"));
		assertFalse(allAll.containsUnusualAA("RT", 69, "_"));
		assertTrue(allAll.containsUnusualAA(Gene.PR, 69, "_"));
		assertFalse(allAll.containsUnusualAA(Gene.RT, 67, "-"));
		assertTrue(allAll.containsUnusualAA("PR", 67, "-"));
	}

}
