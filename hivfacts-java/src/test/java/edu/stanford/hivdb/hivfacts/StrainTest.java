package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import org.junit.Test;

public class StrainTest {

	@Test
	public void testGetDisplayText() {
		assertEquals("HIV-1", Strain.HIV1.getDisplayText());
		assertEquals("HIV-2 Group A", Strain.HIV2A.getDisplayText());
		assertEquals("HIV-2 Group B", Strain.HIV2B.getDisplayText());
	}

	@Test
	public void testGetNucaminoProfile() {
		assertEquals("hiv1b", Strain.HIV1.getNucaminoProfile());
		assertEquals("hiv2a", Strain.HIV2A.getNucaminoProfile());
		assertEquals("hiv2b", Strain.HIV2B.getNucaminoProfile());
	}

}
