package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import org.junit.Test;

public class HIVAPOBECMutationTest {

	@Test
	public void test() {
		HIVAPOBECMutation mut = new HIVAPOBECMutation("RT", 93, 'K');
		assertEquals(HIVGene.valueOf("HIV1RT"), mut.getGene());
		assertEquals(new Integer(93), mut.getPosition());
		assertEquals(new Character('K'), mut.getAA());
		assertEquals(new HIVGenePosition(HIVGene.valueOf("HIV1RT"), 93), mut.getGenePosition());
		assertEquals("HIV1RT:93K", mut.toString());
	}

}
