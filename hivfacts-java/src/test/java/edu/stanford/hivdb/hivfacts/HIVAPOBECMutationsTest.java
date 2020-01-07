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
		HIVAPOBECS inst1 = HIVAPOBECS.getInstance();
		HIVAPOBECS inst2 = HIVAPOBECS.getInstance();
		assertEquals("same singleton instance", inst1, inst2);
	}

	@Test
	public void testloadAPOBECResFailed() {
		expectedEx.expect(ExceptionInInitializerError.class);
		expectedEx.expectMessage("Invalid resource name (apobec/apobec.json)");
		HIVAPOBECS.loadMutationSetFromRes("apobec/apobec.json");
	}

	@Test
	public void testIsApobecMutation() {
		HIVAPOBECS inst = HIVAPOBECS.getInstance();
		assertTrue(inst.isApobecMutation(HIVGene.valueOf("HIV1PR"), 27, 'E'));
		assertFalse(inst.isApobecMutation(HIVGene.valueOf("HIV1PR"), 27, 'G'));
		assertFalse(inst.isApobecMutation(HIVGene.valueOf("HIV1IN"), 263, 'K'));
	}

	@Test
	public void testIsApobecDRM() {
		HIVAPOBECS inst = HIVAPOBECS.getInstance();
		assertFalse(inst.isApobecDRM(HIVGene.valueOf("HIV1PR"), 27, 'E'));
		assertFalse(inst.isApobecDRM(HIVGene.valueOf("HIV1PR"), 27, 'G'));
		assertTrue(inst.isApobecDRM(HIVGene.valueOf("HIV1IN"), 263, 'K'));
		assertFalse(inst.isApobecDRM(HIVGene.valueOf("HIV1IN"), 263, 'R'));
	}

	@Test
	public void testGetApobecMutations() {
		HIVAPOBECS inst = HIVAPOBECS.getInstance();
		assertNotNull(inst.getApobecMutations());
	}

	@Test
	public void testGetApobecDRMs() {
		HIVAPOBECS inst = HIVAPOBECS.getInstance();
		assertNotNull(inst.getApobecDRMs());
	}

}
