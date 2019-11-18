package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class HIVAPOBECMutationsTest {

	@Rule
	public ExpectedException expectedEx = ExpectedException.none();

	@Test
	public void testGetInstanceSuccess() {
		HIVAPOBECMutations inst1 = HIVAPOBECMutations.getInstance();
		HIVAPOBECMutations inst2 = HIVAPOBECMutations.getInstance();
		assertEquals("same singleton instance", inst1, inst2);
	}

	@Test
	public void testloadAPOBECResFailed() {
		expectedEx.expect(ExceptionInInitializerError.class);
		expectedEx.expectMessage("Invalid resource name (apobec/apobec.json)");
		HIVAPOBECMutations.loadAPOBECRes("apobec/apobec.json");
	}

	@Test
	public void testIsApobecMutation() {
		HIVAPOBECMutations inst = HIVAPOBECMutations.getInstance();
		assertTrue(inst.isApobecMutation(HIVGene.valueOf("HIV1PR"), 27, 'E'));
		assertFalse(inst.isApobecMutation(HIVGene.valueOf("HIV1PR"), 27, 'G'));
		assertFalse(inst.isApobecMutation(HIVGene.valueOf("HIV1IN"), 263, 'K'));
	}

	@Test
	public void testIsApobecDRM() {
		HIVAPOBECMutations inst = HIVAPOBECMutations.getInstance();
		assertFalse(inst.isApobecDRM(HIVGene.valueOf("HIV1PR"), 27, 'E'));
		assertFalse(inst.isApobecDRM(HIVGene.valueOf("HIV1PR"), 27, 'G'));
		assertTrue(inst.isApobecDRM(HIVGene.valueOf("HIV1IN"), 263, 'K'));
		assertFalse(inst.isApobecDRM(HIVGene.valueOf("HIV1IN"), 263, 'R'));
	}

	@Test
	public void testGetApobecMutations() {
		HIVAPOBECMutations inst = HIVAPOBECMutations.getInstance();
		assertNotNull(inst.getApobecMutations());
	}

	@Test
	public void testGetApobecDRMs() {
		HIVAPOBECMutations inst = HIVAPOBECMutations.getInstance();
		assertNotNull(inst.getApobecDRMs());
	}

}
