package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import org.junit.Test;

import edu.stanford.hivdb.mutations.AminoAcidPercent;

public class HIVAminoAcidPercentTest {

	@Test
	public void test() {
		AminoAcidPercent mPR5A = new AminoAcidPercent(
			"PR", 5, 'A',
			/* percent   = */ 6.7597712493409226e-06,
			/* count     = */ 1,
			/* total     = */ 147934,
			/* reason    = */ "PCNT",
			/* isUnusual = */ true
		);
		assertEquals(HIVGene.valueOf("HIV1PR"), mPR5A.getGene());
		assertEquals(5, (int) mPR5A.position);
		assertEquals('A', (char) mPR5A.aa);
		assertEquals(6.7597712493409226e-06, mPR5A.percent, 1e-18);
		assertEquals(1, (int) mPR5A.count);
		assertEquals(147934, (int) mPR5A.total);
		assertEquals("PCNT", mPR5A.reason);
		assertTrue(mPR5A.isUnusual);
	}

}
