package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import org.junit.Test;

public class MutTypeTest {

	@Test
	public void testIsDRMType() {
		assertTrue(MutType.Major.isDRMType());
		assertTrue(MutType.NRTI.isDRMType());
		assertTrue(MutType.NNRTI.isDRMType());
		assertFalse(MutType.Accessory.isDRMType());
		assertFalse(MutType.Other.isDRMType());
	}

}
