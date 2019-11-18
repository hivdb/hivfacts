package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;
import java.util.HashSet;
import java.util.Set;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import com.google.common.collect.Sets;

public class GenePositionTest {
	static final int MAX_PR_POS = 99;
	static final int MAX_RT_POS = 560;
	static final int MAX_IN_POS = 288;
	static final int TOTAL_POS = MAX_PR_POS + MAX_RT_POS + MAX_IN_POS;

	@Rule
	public ExpectedException expectedEx = ExpectedException.none();


	// Potential issues:
	// 1. Doesn't throw out of bounds exceptions.
	@Test
	public void testConstruction() {
		final HIVGenePosition prGP = new HIVGenePosition(HIVGene.valueOf("HIV1PR"), MAX_PR_POS);
		final HIVGenePosition rtGP = new HIVGenePosition(HIVGene.valueOf("HIV1RT"), MAX_RT_POS);
		final HIVGenePosition inGP = new HIVGenePosition(HIVGene.valueOf("HIV1IN"), MAX_IN_POS);
		final HIVGenePosition inGPMin = new HIVGenePosition(HIVGene.valueOf("HIV1IN"), Integer.MIN_VALUE);
		final HIVGenePosition inGPMax = new HIVGenePosition(HIVGene.valueOf("HIV1IN"), Integer.MAX_VALUE);
		assertEquals(prGP.gene, HIVGene.valueOf("HIV1PR"));
		assertEquals(rtGP.gene, HIVGene.valueOf("HIV1RT"));
		assertEquals(inGP.gene, HIVGene.valueOf("HIV1IN"));
		assertEquals(prGP.position, Integer.valueOf(MAX_PR_POS));
		assertEquals(rtGP.position, Integer.valueOf(MAX_RT_POS));
		assertEquals(inGP.position, Integer.valueOf(MAX_IN_POS));
		assertEquals(inGPMax.position, Integer.valueOf(Integer.MAX_VALUE));
		assertEquals(inGPMin.position, Integer.valueOf(Integer.MIN_VALUE));
	}

	@Test
	public void testGetGenePositionsBetween() {
		assertEquals(
			Sets.newHashSet(
				new HIVGenePosition("HIV1PR:97"),
				new HIVGenePosition("HIV1PR:98"),
				new HIVGenePosition("HIV1PR:99"),
				new HIVGenePosition("HIV1RT:1"),
				new HIVGenePosition("HIV1RT:2")
			), HIVGenePosition.getGenePositionsBetween(
				new HIVGenePosition("HIV1PR:97"),
				new HIVGenePosition("HIV1RT:2")
			));

		Set<HIVGenePosition> gps1 = HIVGenePosition.getGenePositionsBetween(
			new HIVGenePosition("HIV1PR:50"),
			new HIVGenePosition("HIV1IN:50")
		);
		assertEquals(660, gps1.size());
		Set<HIVGenePosition> gps2 = HIVGenePosition.getGenePositionsBetween(
			new HIVGenePosition("HIV2APR:50"),
			new HIVGenePosition("HIV2AIN:50")
		);
		assertEquals(659, gps2.size());
		Set<HIVGenePosition> gps3 = HIVGenePosition.getGenePositionsBetween(
			new HIVGenePosition("HIV1PR:50"),
			new HIVGenePosition("HIV1PR:95")
		);
		assertEquals(46, gps3.size());
		Set<HIVGenePosition> gps4 = HIVGenePosition.getGenePositionsBetween(
			new HIVGenePosition("HIV2BPR:50"),
			new HIVGenePosition("HIV2BRT:50")
		);
		assertEquals(100, gps4.size());
	}
	
	@Test
	public void testGetPolPosition() {
		assertEquals(1, (int) new HIVGenePosition("HIV1PR:1").getPolPosition());
		assertEquals(100, (int) new HIVGenePosition("HIV1RT:1").getPolPosition());
		assertEquals(660, (int) new HIVGenePosition("HIV1IN:1").getPolPosition());
		assertEquals(444, (int) new HIVGenePosition("HIV2ART:345").getPolPosition());
		assertEquals(445, (int) new HIVGenePosition("HIV2ART:346").getPolPosition());
		assertEquals(658, (int) new HIVGenePosition("HIV2ART:559").getPolPosition());
		assertEquals(659, (int) new HIVGenePosition("HIV2AIN:1").getPolPosition());
	}

	@Test
	public void testGetGenePositionsBetweenMismatchedStrain() {
		expectedEx.expect(IllegalArgumentException.class);
		expectedEx.expectMessage("Virus strain of `start` and `end` positions must be the same.");
		HIVGenePosition.getGenePositionsBetween(
			new HIVGenePosition("HIV1PR:50"),
			new HIVGenePosition("HIV2AIN:50")
		);
	}

	// Potential issues:
	// 1. Case sensitive
	// 2. Doesn't throw exceptions for malformed strings.
	//	  Perhaps this method could parse the string with a regex and
	//	  throw an exception if the input doesn't match it.
	@Test
	public void testConstructionFromString() {
		final HIVGenePosition prGP = new HIVGenePosition("HIV1PR:99");
		final HIVGenePosition rtGP = new HIVGenePosition("HIV1RT:560");
		final HIVGenePosition inGP = new HIVGenePosition("HIV1IN:288");
		assertEquals(prGP.gene, HIVGene.valueOf("HIV1PR"));
		assertEquals(rtGP.gene, HIVGene.valueOf("HIV1RT"));
		assertEquals(inGP.gene, HIVGene.valueOf("HIV1IN"));
		assertEquals(prGP.position, Integer.valueOf(MAX_PR_POS));
		assertEquals(rtGP.position, Integer.valueOf(MAX_RT_POS));
		assertEquals(inGP.position, Integer.valueOf(MAX_IN_POS));
		assertEquals("HIV1PR:99", prGP.toString());
		assertEquals("HIV1RT:560", rtGP.toString());
		assertEquals("HIV1IN:288", inGP.toString());
	}

	@Test(expected = NullPointerException.class)
	public void testCompareToException() {
		final HIVGenePosition gp = new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 1);
		assertEquals(gp.compareTo(null), 1);
	}

	@Test
	@SuppressWarnings("unlikely-arg-type")
	public void testEquals() {
		HIVGenePosition gpPR5 = new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 5);
		assertFalse(gpPR5.equals(null));
		assertTrue(gpPR5.equals(gpPR5));
		assertFalse(gpPR5.equals("HIV1PR:5"));
		assertTrue(gpPR5.equals(new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 5)));
		assertFalse(gpPR5.equals(new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 6)));
		assertFalse(gpPR5.equals(new HIVGenePosition(HIVGene.valueOf("HIV1RT"), 5)));
		assertFalse(gpPR5.equals(new HIVGenePosition(HIVGene.valueOf("HIV1RT"), 6)));

		final HIVGenePosition prGP = new HIVGenePosition(HIVGene.valueOf("HIV1PR"), MAX_PR_POS);
		final HIVGenePosition prGPFromStr = new HIVGenePosition("HIV1PR:99");
		final HIVGenePosition rtGPFromStr = new HIVGenePosition("HIV1RT:560");
		assertTrue(prGPFromStr.equals(prGPFromStr));
		assertTrue(prGPFromStr.equals(prGP));
		assertFalse(prGPFromStr.equals(null));
		assertFalse(prGPFromStr.equals(rtGPFromStr));
		assertFalse(prGPFromStr.equals(HIVGene.valueOf("HIV1PR")));

	}

	@Test
	public void testCompareTo() {
		assertEquals(0, new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 5).compareTo(new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 5)));
		assertEquals(1, new HIVGenePosition(HIVGene.valueOf("HIV1RT"), 5).compareTo(new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 5)));
		assertEquals(1, new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 6).compareTo(new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 5)));
		expectedEx.expect(NullPointerException.class);
		expectedEx.expectMessage("Null is incomprable.");
		new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 5).compareTo(null);

		final HIVGenePosition prGPMin = new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 1);
		final HIVGenePosition prGPMid = new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 50);
		final HIVGenePosition prGPMax = new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 99);
		final HIVGenePosition rtGP = new HIVGenePosition(HIVGene.valueOf("HIV1RT"), 99);
		assertEquals(prGPMin.compareTo(rtGP), -1);
		assertEquals(prGPMin.compareTo(prGPMin), 0);
		assertEquals(prGPMin.compareTo(prGPMid), -1);
		assertEquals(prGPMax.compareTo(prGPMid), 1);
	}

	@Test
	public void testToString() {
		assertEquals("HIV1PR:5", new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 5).toString());
	}

	@Test
	public void testHashCode() {
		assertEquals(1, new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 6).hashCode() - new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 5).hashCode());
		assertEquals(new HIVGenePosition(HIVGene.valueOf("HIV1RT"), 5).hashCode(), new HIVGenePosition(HIVGene.valueOf("HIV1RT"), 5).hashCode());
		assertNotEquals(new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 5).hashCode(), new HIVGenePosition(HIVGene.valueOf("HIV1RT"), 5).hashCode());

		Set<Integer> hashCodes = new HashSet<Integer>();
		for (int pos = 1; pos <= MAX_RT_POS; pos++) {
			if (pos <= MAX_PR_POS) hashCodes.add(new HIVGenePosition(HIVGene.valueOf("HIV1PR"), pos).hashCode());
			if (pos <= MAX_IN_POS) hashCodes.add(new HIVGenePosition(HIVGene.valueOf("HIV1IN"), pos).hashCode());
			hashCodes.add(new HIVGenePosition(HIVGene.valueOf("HIV1RT"), pos).hashCode());
		}
		assertEquals(hashCodes.size(), TOTAL_POS);
	}
}
