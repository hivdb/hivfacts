package edu.stanford.hivdb.aapcnt;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class HIVAminoAcidPercentsTest {
	
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

}