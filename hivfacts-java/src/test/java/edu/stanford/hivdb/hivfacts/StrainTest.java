package edu.stanford.hivdb.hivfacts;

import static org.junit.Assert.*;

import org.junit.Test;

public class StrainTest {

	@Test
	public void testGetDisplayText() {
		assertEquals("HIV-1", HIVStrain.HIV1.getDisplayText());
		assertEquals("HIV-2 Group A", HIVStrain.HIV2A.getDisplayText());
		assertEquals("HIV-2 Group B", HIVStrain.HIV2B.getDisplayText());
	}

	@Test
	public void testGetNucaminoProfile() {
		assertEquals("hiv1b", HIVStrain.HIV1.getNucaminoProfile());
		assertEquals("hiv2a", HIVStrain.HIV2A.getNucaminoProfile());
		assertEquals("hiv2b", HIVStrain.HIV2B.getNucaminoProfile());
	}

}
