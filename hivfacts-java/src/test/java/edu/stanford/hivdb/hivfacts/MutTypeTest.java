package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import org.junit.Test;

import edu.stanford.hivdb.mutations.MutationType;

public class MutTypeTest {

	@Test
	public void testIsDRMType() {
		assertTrue(MutationType.Major.isDRMType());
		assertTrue(MutationType.NRTI.isDRMType());
		assertTrue(MutationType.NNRTI.isDRMType());
		assertFalse(MutationType.Accessory.isDRMType());
		assertFalse(MutationType.Other.isDRMType());
	}

}
