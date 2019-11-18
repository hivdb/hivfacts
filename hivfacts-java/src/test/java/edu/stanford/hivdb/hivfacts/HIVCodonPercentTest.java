package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import org.junit.Test;

public class HIVCodonPercentTest {

	@Test
	public void test() {
		HIVCodonPercent mPR5A = new HIVCodonPercent(
			"PR", 5, "GCT", 'A',
			/* percent   = */ 5.664567000498482e-06,
			/* count     = */ 1,
			/* total     = */ 176536
		);
		assertEquals(HIVGene.valueOf("HIV1PR"), mPR5A.getGene());
		assertEquals(5, (int) mPR5A.position);
		assertEquals("GCT", (String) mPR5A.codon);
		assertEquals('A', (char) mPR5A.aa);
		assertEquals(5.664567000498482e-06, mPR5A.percent, 1e-18);
		assertEquals(1, (int) mPR5A.count);
		assertEquals(176536, (int) mPR5A.total);
	}

}
