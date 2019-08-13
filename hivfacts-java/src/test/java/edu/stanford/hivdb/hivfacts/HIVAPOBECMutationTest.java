package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import org.junit.Test;

public class HIVAPOBECMutationTest {

	@Test
	public void test() {
		HIVAPOBECMutation mut = new HIVAPOBECMutation("RT", 93, 'K');
		assertEquals("RT", mut.getGene());
		assertEquals(new Integer(93), mut.getPosition());
		assertEquals(new Character('K'), mut.getAA());
		assertEquals(new GenePosition("RT", 93), mut.getGenePosition());
		assertEquals("RT:93K", mut.toString());
	}

}
