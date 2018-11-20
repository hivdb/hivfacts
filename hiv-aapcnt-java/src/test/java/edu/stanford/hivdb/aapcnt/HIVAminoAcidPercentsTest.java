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
		HIVAminoAcidPercents naiveAll01 = HIVAminoAcidPercents.getInstance("naive", "All");
		HIVAminoAcidPercents naiveAll02 = HIVAminoAcidPercents.getInstance("naive", "All");
		assertEquals("same singleton instance", naiveAll01, naiveAll02);
	}
	
	@Test
	public void testGetInstanceFailCase1() {
		expectedEx.expect(ExceptionInInitializerError.class);
		expectedEx.expectMessage("Invalid resource name (naiveE.json)");
		HIVAminoAcidPercents.getInstance("naive", "E");
	}

	@Test
	public void testGetInstanceFailCase2() {
		expectedEx.expect(ExceptionInInitializerError.class);
		expectedEx.expectMessage("Invalid resource name (aaaaaaaAll.json)");
		HIVAminoAcidPercents.getInstance("aaaaaaa", "All");
	}
	
	@Test
	public void testGsonLoad() {
		HIVAminoAcidPercents naiveAll = HIVAminoAcidPercents.getInstance("naive", "All");
		assertEquals((99 + 560 + 288) * 23, naiveAll.aminoAcidPcnts.size());
		HIVAminoAcidPercent mutPR1A = naiveAll.aminoAcidPcnts.get(0);
		assertEquals("PR", mutPR1A.gene);
		assertEquals(1, (int) mutPR1A.position);
		assertEquals('A', (char) mutPR1A.aa);
	}
	
	@Test
	public void testGetAll() {
		HIVAminoAcidPercents naiveAll = HIVAminoAcidPercents.getInstance("naive", "All");
		List<HIVAminoAcidPercent> allAAPcnts = naiveAll.get();
		assertEquals(naiveAll.aminoAcidPcnts, allAAPcnts); // equal values
		assertFalse(naiveAll.aminoAcidPcnts == allAAPcnts); // different references
		assertEquals((99 + 560 + 288) * 23, allAAPcnts.size());
	}
	
	@Test
	public void testGetByGene() {
		HIVAminoAcidPercents naiveAll = HIVAminoAcidPercents.getInstance("naive", "All");
		List<HIVAminoAcidPercent> gRTAAPcnts = naiveAll.get("RT");
		List<HIVAminoAcidPercent> gRTAAPcntsByEnum = naiveAll.get(Gene.RT);
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
		HIVAminoAcidPercents naiveAll = HIVAminoAcidPercents.getInstance("naive", "All");
		List<HIVAminoAcidPercent> gpIN263AAPcnts = naiveAll.get("IN", 263);
		List<HIVAminoAcidPercent> gpIN263AAPcntsByEnum = naiveAll.get(Gene.IN, 263);
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
		HIVAminoAcidPercents naiveAll = HIVAminoAcidPercents.getInstance("naive", "All");
		HIVAminoAcidPercent mutIN263R = naiveAll.get("IN", 263, 'R');
		HIVAminoAcidPercent mutIN263RByEnum = naiveAll.get(Gene.IN, 263, 'R');
		assertEquals(mutIN263R, mutIN263RByEnum); // equal value
		assertTrue(mutIN263R == mutIN263RByEnum); // same references
		assertEquals("IN", mutIN263R.gene);
		assertEquals(263, (int) mutIN263R.position);
		assertEquals('R', (char) mutIN263R.aa);
	}
	
	@Test
	public void testGetHighestPercentValue() {
		HIVAminoAcidPercents naiveAll = HIVAminoAcidPercents.getInstance("naive", "All");
		double highestVal = naiveAll.getHighestAAPercentValue("PR", 23, "AHI");
		double highestValByEnum = naiveAll.getHighestAAPercentValue(Gene.PR, 23, "AHI");
		assertEquals(highestVal, highestValByEnum, 1e-18);
		double expectedHighestVal = .0;
		for (char aa : "AHI".toCharArray()) {
			double pcntVal = naiveAll.get("PR", 23, aa).percent;
			expectedHighestVal = Math.max(expectedHighestVal, pcntVal);
		}
		assertEquals(expectedHighestVal, highestVal, 1e-18);
		// This is intended to fail for every version update.
		// You must manually check and correct this number.
		assertEquals(0.002121155805622665, highestVal, 1e-18);
	}
	
	@Test
	public void testContainsUnusualAA() {
		HIVAminoAcidPercents naiveAll = HIVAminoAcidPercents.getInstance("naive", "All");
		assertTrue(naiveAll.containsUnusualAA("RT", 5, "I*"));
		assertFalse(naiveAll.containsUnusualAA(Gene.PR, 6, "W*")); // APOBEC stop codon
	}

}