package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class HIVCodonPercentsTest {

	@Rule
	public ExpectedException expectedEx = ExpectedException.none();

	@Test
	public void testGetInstanceSuccess() {
		HIVCodonPercents allall01 = HIVCodonPercents.getInstance("all", "all");
		HIVCodonPercents allall02 = HIVCodonPercents.getInstance("all", "all");
		assertEquals("same singleton instance", allall01, allall02);
	}

	@Test
	public void testGetInstanceFailCase1() {
		expectedEx.expect(ExceptionInInitializerError.class);
		expectedEx.expectMessage("Invalid resource name (codonpcnt/rx-all_subtype-E.json)");
		HIVCodonPercents.getInstance("all", "E");
	}

	@Test
	public void testGetInstanceFailCase2() {
		expectedEx.expect(ExceptionInInitializerError.class);
		expectedEx.expectMessage("Invalid resource name (codonpcnt/rx-aaaaaaa_subtype-all.json)");
		HIVCodonPercents.getInstance("aaaaaaa", "all");
	}

	@Test
	public void testGsonLoad() {
		HIVCodonPercents allall = HIVCodonPercents.getInstance("all", "all");
		assertNotNull(allall.codonPcnts);
		HIVCodonPercent mutPR1ACA = allall.codonPcnts.get(0);
		assertEquals(Gene.valueOf("HIV1PR"), mutPR1ACA.getGene());
		assertEquals(new Integer(1), mutPR1ACA.position);
		assertEquals(new Character('T'), mutPR1ACA.aa);
		assertEquals("ACA", mutPR1ACA.codon);
	}

	@Test
	public void testGetall() {
		HIVCodonPercents allall = HIVCodonPercents.getInstance("all", "all");
		List<HIVCodonPercent> allCdPcnts = allall.get();
		assertEquals(allall.codonPcnts, allCdPcnts); // equal values
		assertFalse(allall.codonPcnts == allCdPcnts); // different references
	}

	@Test
	public void testGetByGene() {
		HIVCodonPercents allall = HIVCodonPercents.getInstance("all", "all");
		List<HIVCodonPercent> gRTCdPcnts = allall.get(Gene.valueOf("HIV1RT"));
		int i = 0;
		int pos = 1;
		char[] acgt = "ACGT".toCharArray();
		// Some of the codons might never occur. However the list
		// should be still in alphabetical order.
		for (char na1 : acgt) {
			for (char na2 : acgt) {
				for (char na3 : acgt) {
					if (i >= gRTCdPcnts.size()) {
						return;
					}
					HIVCodonPercent mut = gRTCdPcnts.get(i ++);
					assertEquals(Gene.valueOf("HIV1RT"), mut.getGene());
					assertEquals(pos, (int) mut.position);
					String expectedCodon = "" + na1 + na2 + na3;
					assertTrue(mut.codon.compareTo(expectedCodon) >= 0);
					if (mut.codon.compareTo(expectedCodon) != 0) {
						i --;
					}
					if (expectedCodon.equals("TTT")) {
						pos ++;
					}
				}
			}
		}
	}

	@Test
	public void testGetByGenePos() {
		HIVCodonPercents allall = HIVCodonPercents.getInstance("all", "all");
		List<HIVCodonPercent> gpIN263CdPcnts = allall.get(Gene.valueOf("HIV1IN"), 263);
		int i = 0;
		char[] acgt = "ACGT".toCharArray();
		// Some of the codons might never occur. However the list
		// should be still in alphabetical order.
		for (char na1 : acgt) {
			for (char na2 : acgt) {
				for (char na3 : acgt) {
					if (i >= gpIN263CdPcnts.size()) {
						return;
					}
					HIVCodonPercent mut = gpIN263CdPcnts.get(i ++);
					assertEquals(Gene.valueOf("HIV1IN"), mut.getGene());
					assertEquals(263, (int) mut.position);
					assertTrue(mut.codon.compareTo("" + na1 + na2 + na3) >= 0);
					if (mut.codon.compareTo("" + na1 + na2 + na3) != 0) {
						i --;
					}
				}
			}
		}
	}

	@Test
	public void testGetByMut() {
		HIVCodonPercents allall = HIVCodonPercents.getInstance("all", "all");
		HIVCodonPercent mutIN263AGG = allall.get(Gene.valueOf("HIV1IN"), 263, "AGG");
		assertEquals(Gene.valueOf("HIV1IN"), mutIN263AGG.getGene());
		assertEquals(263, (int) mutIN263AGG.position);
		assertEquals('R', (char) mutIN263AGG.aa);
		assertEquals("AGG", mutIN263AGG.codon);
		HIVCodonPercent mutIN5ins = allall.get(Gene.valueOf("HIV1IN"), 5, "ins");
		assertEquals(Gene.valueOf("HIV1IN"), mutIN5ins.getGene());
		assertEquals(5, (int) mutIN5ins.position);
		assertEquals('_', (int) mutIN5ins.aa);
		assertEquals("ins", mutIN5ins.codon);
		HIVCodonPercent mutIN5del = allall.get(Gene.valueOf("HIV1IN"), 5, "del");
		assertEquals(Gene.valueOf("HIV1IN"), mutIN5del.getGene());
		assertEquals(5, (int) mutIN5del.position);
		assertEquals('-', (int) mutIN5del.aa);
		assertEquals("del", mutIN5del.codon);
	}

	@Test
	public void testGetByMutNotOccur() {
		HIVCodonPercents allall = HIVCodonPercents.getInstance("all", "all");
		HIVCodonPercent mutIN263GTG = allall.get(Gene.valueOf("HIV1IN"), 263, "GTG");
		assertEquals(Gene.valueOf("HIV1IN"), mutIN263GTG.getGene());
		assertEquals(263, (int) mutIN263GTG.position);
		// TODO: don't know the translation. Probably should introduce codon translator
		assertEquals('X', (char) mutIN263GTG.aa);
		assertEquals("GTG", mutIN263GTG.codon);
	}

	@Test
	public void testGetByInvalidMut() {
		HIVCodonPercents allall = HIVCodonPercents.getInstance("all", "all");
		expectedEx.expect(IllegalArgumentException.class);
		expectedEx.expectMessage("Invalid argument codon \"EEE\" at HIV1IN1");
		allall.get(Gene.valueOf("HIV1IN"), 1, "EEE");
	}
	
	@Test
	public void testGetByInvalidMut2() {
		HIVCodonPercents allall = HIVCodonPercents.getInstance("all", "all");
		expectedEx.expect(IllegalArgumentException.class);
		expectedEx.expectMessage("Invalid argument codon \"\" at HIV1IN1");
		allall.get(Gene.valueOf("HIV1IN"), 1, "");
	}

	@Test
	public void testGetOutOfRange() {
		HIVCodonPercents allall = HIVCodonPercents.getInstance("all", "all");
		expectedEx.expect(IllegalArgumentException.class);
		expectedEx.expectMessage("Argument 'pos' is out of range: 100");
		allall.get(Gene.valueOf("HIV1PR"), 100, "AAA");

	}

	@Test
	public void testGetHighestPercentValue() {
		HIVCodonPercents allall = HIVCodonPercents.getInstance("all", "all");
		double highestVal = allall.getHighestCodonPercentValue(Gene.valueOf("HIV1PR"), 23, "AAA", "AAG", "ATA", "ATC");
		double expectedHighestVal = .0;
		for (String codon : new String[] {"AAA", "AAG", "ATA", "ATC"}) {
			double pcntVal = allall.get(Gene.valueOf("HIV1PR"), 23, codon).percent;
			expectedHighestVal = Math.max(expectedHighestVal, pcntVal);
		}
		assertEquals(expectedHighestVal, highestVal, 1e-18);

		// These are intended to fail for every version update.
		// You must manually check and correct these numbers.
		assertEquals(0.00166933, highestVal, 1e-8);
		assertEquals(0.76700478, allall.getHighestCodonPercentValue(Gene.valueOf("HIV1RT"), 67, "GAC"), 1e-8);
		assertEquals(0.00381235, allall.getHighestCodonPercentValue(Gene.valueOf("HIV1RT"), 69, "AAA", "AGC", "AGT"), 1e-8);
		assertEquals(0.0, allall.getHighestCodonPercentValue(Gene.valueOf("HIV1RT"), 67, "TGG"), 1e-8);
		assertEquals(0.06986288, allall.getHighestCodonPercentValue(Gene.valueOf("HIV1RT"), 67, "AAC", "TGA"), 1e-8);
		assertEquals(0.0, allall.getHighestCodonPercentValue(Gene.valueOf("HIV1RT"), 67, "TGA"), 1e-8);

		expectedEx.expect(IllegalArgumentException.class);
		expectedEx.expectMessage("Invalid argument codon \"EEE\" at HIV1IN1");
		allall.getHighestCodonPercentValue(Gene.valueOf("HIV1IN"), 1, "EEE");
	}

}
