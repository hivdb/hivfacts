package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class GenePositionTest {

	@Rule
	public ExpectedException expectedEx = ExpectedException.none();

	@Test
	@SuppressWarnings("unlikely-arg-type")
	public void testEquals() {
		GenePosition gpPR5 = new GenePosition("PR", 5);
		assertFalse(gpPR5.equals(null));
		assertTrue(gpPR5.equals(gpPR5));
		assertFalse(gpPR5.equals("PR:5"));
		assertTrue(gpPR5.equals(new GenePosition("PR", 5)));
		assertFalse(gpPR5.equals(new GenePosition("PR", 6)));
		assertFalse(gpPR5.equals(new GenePosition("RT", 5)));
		assertFalse(gpPR5.equals(new GenePosition("RT", 6)));
		
	}
	
	@Test
	public void testCompareTo() {
		assertEquals(0, new GenePosition("PR", 5).compareTo(new GenePosition("PR", 5)));
		assertEquals(2, new GenePosition("RT", 5).compareTo(new GenePosition("PR", 5)));
		assertEquals(1, new GenePosition("PR", 6).compareTo(new GenePosition("PR", 5)));
		expectedEx.expect(NullPointerException.class);
		expectedEx.expectMessage("Null is incomprable.");
		new GenePosition("PR", 5).compareTo(null);
	}
	
	@Test
	public void testToString() {
		assertEquals("PR:5", new GenePosition("PR", 5).toString());
	}

	@Test
	public void testHashCode() {
		assertEquals(1708231282, new GenePosition("PR", 5).hashCode());
		assertEquals(1708231283, new GenePosition("PR", 6).hashCode());
		assertEquals(-1594273870, new GenePosition("RT", 5).hashCode());
	}

}
