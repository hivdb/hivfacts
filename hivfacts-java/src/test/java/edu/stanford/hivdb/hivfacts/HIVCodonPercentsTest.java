package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import edu.stanford.hivdb.mutations.CodonPercent;
import edu.stanford.hivdb.mutations.CodonPercents;

public class HIVCodonPercentsTest {

	@Rule
	public ExpectedException expectedEx = ExpectedException.none();

	@Test
	public void testGetInstanceSuccess() {
		CodonPercents allall01 = CodonPercents.getInstance("all", "all");
		CodonPercents allall02 = CodonPercents.getInstance("all", "all");
		assertEquals("same singleton instance", allall01, allall02);
	}

	@Test
	public void testGetInstanceFailCase1() {
		expectedEx.expect(ExceptionInInitializerError.class);
		expectedEx.expectMessage("Invalid resource name (codonpcnt/rx-all_subtype-E.json)");
		CodonPercents.getInstance("all", "E");
	}

	@Test
	public void testGetInstanceFailCase2() {
		expectedEx.expect(ExceptionInInitializerError.class);
		expectedEx.expectMessage("Invalid resource name (codonpcnt/rx-aaaaaaa_subtype-all.json)");
		CodonPercents.getInstance("aaaaaaa", "all");
	}

	@Test
	public void testGsonLoad() {
		CodonPercents allall = CodonPercents.getInstance("all", "all");
		assertNotNull(allall.codonPcnts);
		CodonPercent mutPR1ACA = allall.codonPcnts.get(0);
		assertEquals(HIVGene.valueOf("HIV1PR"), mutPR1ACA.getGene());
		assertEquals(new Integer(1), mutPR1ACA.position);
		assertEquals(new Character('T'), mutPR1ACA.aa);
		assertEquals("ACA", mutPR1ACA.codon);
	}

	@Test
	public void testGetall() {
		CodonPercents allall = CodonPercents.getInstance("all", "all");
		List<CodonPercent> allCdPcnts = allall.get();
		assertEquals(allall.codonPcnts, allCdPcnts); // equal values
		assertFalse(allall.codonPcnts == allCdPcnts); // different references
	}

	@Test
	public void testGetByGene() {
		CodonPercents allall = CodonPercents.getInstance("all", "all");
		List<CodonPercent> gRTCdPcnts = allall.get(HIVGene.valueOf("HIV1RT"));
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
					CodonPercent mut = gRTCdPcnts.get(i ++);
					assertEquals(HIVGene.valueOf("HIV1RT"), mut.getGene());
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
		CodonPercents allall = CodonPercents.getInstance("all", "all");
		List<CodonPercent> gpIN263CdPcnts = allall.get(HIVGene.valueOf("HIV1IN"), 263);
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
					CodonPercent mut = gpIN263CdPcnts.get(i ++);
					assertEquals(HIVGene.valueOf("HIV1IN"), mut.getGene());
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
		CodonPercents allall = CodonPercents.getInstance("all", "all");
		CodonPercent mutIN263AGG = allall.get(HIVGene.valueOf("HIV1IN"), 263, "AGG");
		assertEquals(HIVGene.valueOf("HIV1IN"), mutIN263AGG.getGene());
		assertEquals(263, (int) mutIN263AGG.position);
		assertEquals('R', (char) mutIN263AGG.aa);
		assertEquals("AGG", mutIN263AGG.codon);
		CodonPercent mutIN5ins = allall.get(HIVGene.valueOf("HIV1IN"), 5, "ins");
		assertEquals(HIVGene.valueOf("HIV1IN"), mutIN5ins.getGene());
		assertEquals(5, (int) mutIN5ins.position);
		assertEquals('_', (int) mutIN5ins.aa);
		assertEquals("ins", mutIN5ins.codon);
		CodonPercent mutIN5del = allall.get(HIVGene.valueOf("HIV1IN"), 5, "del");
		assertEquals(HIVGene.valueOf("HIV1IN"), mutIN5del.getGene());
		assertEquals(5, (int) mutIN5del.position);
		assertEquals('-', (int) mutIN5del.aa);
		assertEquals("del", mutIN5del.codon);
	}

	@Test
	public void testGetByMutNotOccur() {
		CodonPercents allall = CodonPercents.getInstance("all", "all");
		CodonPercent mutIN263GTG = allall.get(HIVGene.valueOf("HIV1IN"), 263, "GTG");
		assertEquals(HIVGene.valueOf("HIV1IN"), mutIN263GTG.getGene());
		assertEquals(263, (int) mutIN263GTG.position);
		// TODO: don't know the translation. Probably should introduce codon translator
		assertEquals('X', (char) mutIN263GTG.aa);
		assertEquals("GTG", mutIN263GTG.codon);
	}

	@Test
	public void testGetByInvalidMut() {
		CodonPercents allall = CodonPercents.getInstance("all", "all");
		expectedEx.expect(IllegalArgumentException.class);
		expectedEx.expectMessage("Invalid argument codon \"EEE\" at HIV1IN1");
		allall.get(HIVGene.valueOf("HIV1IN"), 1, "EEE");
	}
	
	@Test
	public void testGetByInvalidMut2() {
		CodonPercents allall = CodonPercents.getInstance("all", "all");
		expectedEx.expect(IllegalArgumentException.class);
		expectedEx.expectMessage("Invalid argument codon \"\" at HIV1IN1");
		allall.get(HIVGene.valueOf("HIV1IN"), 1, "");
	}

	@Test
	public void testGetOutOfRange() {
		CodonPercents allall = CodonPercents.getInstance("all", "all");
		expectedEx.expect(IllegalArgumentException.class);
		expectedEx.expectMessage("Argument 'pos' is out of range: 100");
		allall.get(HIVGene.valueOf("HIV1PR"), 100, "AAA");

	}

	@Test
	public void testGetHighestPercentValue() {
		CodonPercents allall = CodonPercents.getInstance("all", "all");
		double highestVal = allall.getHighestCodonPercentValue(HIVGene.valueOf("HIV1PR"), 23, "AAA", "AAG", "ATA", "ATC");
		double expectedHighestVal = .0;
		for (String codon : new String[] {"AAA", "AAG", "ATA", "ATC"}) {
			double pcntVal = allall.get(HIVGene.valueOf("HIV1PR"), 23, codon).percent;
			expectedHighestVal = Math.max(expectedHighestVal, pcntVal);
		}
		assertEquals(expectedHighestVal, highestVal, 1e-18);

		// These are intended to fail for every version update.
		// You must manually check and correct these numbers.
		assertEquals(0.00166933, highestVal, 1e-8);
		assertEquals(0.76700478, allall.getHighestCodonPercentValue(HIVGene.valueOf("HIV1RT"), 67, "GAC"), 1e-8);
		assertEquals(0.00381235, allall.getHighestCodonPercentValue(HIVGene.valueOf("HIV1RT"), 69, "AAA", "AGC", "AGT"), 1e-8);
		assertEquals(0.0, allall.getHighestCodonPercentValue(HIVGene.valueOf("HIV1RT"), 67, "TGG"), 1e-8);
		assertEquals(0.06986288, allall.getHighestCodonPercentValue(HIVGene.valueOf("HIV1RT"), 67, "AAC", "TGA"), 1e-8);
		assertEquals(0.0, allall.getHighestCodonPercentValue(HIVGene.valueOf("HIV1RT"), 67, "TGA"), 1e-8);

		expectedEx.expect(IllegalArgumentException.class);
		expectedEx.expectMessage("Invalid argument codon \"EEE\" at HIV1IN1");
		allall.getHighestCodonPercentValue(HIVGene.valueOf("HIV1IN"), 1, "EEE");
	}

}
