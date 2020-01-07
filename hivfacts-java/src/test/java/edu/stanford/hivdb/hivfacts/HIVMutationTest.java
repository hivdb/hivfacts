/*

    Copyright (C) 2017 Stanford HIVDB team

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.hivfacts;

import org.junit.Test;
import com.google.common.collect.Sets;

import edu.stanford.hivdb.hivfacts.HIVGene;
import edu.stanford.hivdb.hivfacts.HIVGenePosition;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.hivfacts.HIVAAMutation;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class HIVMutationTest {

	@Test
	public void testNormalizeAAChars() {
		assertEquals(null, HIVAAMutation.normalizeAAChars(null));
		assertEquals(Sets.newHashSet('_'), HIVAAMutation.normalizeAAChars(Sets.newHashSet('#')));
		assertEquals(Sets.newHashSet('_'), HIVAAMutation.normalizeAAChars(Sets.newHashSet('i')));
		assertEquals(Sets.newHashSet('-'), HIVAAMutation.normalizeAAChars(Sets.newHashSet('~')));
		assertEquals(Sets.newHashSet('-'), HIVAAMutation.normalizeAAChars(Sets.newHashSet('d')));
		assertEquals(Sets.newHashSet('*'), HIVAAMutation.normalizeAAChars(Sets.newHashSet('Z')));
		assertEquals(Sets.newHashSet('*'), HIVAAMutation.normalizeAAChars(Sets.newHashSet('.')));
	}

	@Test(expected=IllegalArgumentException.class)
	public void testPositionOutOfGene() {
		new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 100, 'A');
	}

	@Test(expected=IllegalArgumentException.class)
	public void testMergesWithNotSameGene() {
		new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, 'A')
			.mergesWith(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 68, 'C'));
	}

	@Test(expected=IllegalArgumentException.class)
	public void testMergesWithNotSamePos() {
		new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, 'A')
			.mergesWith(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 68, 'C'));
	}

	public void testMergesWithIndel() {
		assertEquals(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, new char[] {'A', '-'}),
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, 'A')
				.mergesWith(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, '-')));
	}

	public void testIndelMergesWith() {
		assertEquals(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, new char[] {'A', '-'}),
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, '-')
				.mergesWith(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, 'A')));
	}

	@Test
	public void testMergesWith() {
		assertEquals(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, "AC".toCharArray()),
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, 'A')
				.mergesWith(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, 'C')));
		assertEquals(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, "AC".toCharArray()),
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, "AC".toCharArray())
				.mergesWith(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, 'C')));
	}

	@Test
	public void testSubtractsBy() {
		HIVAAMutation pr67ANXDMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 67, "ANXD".toCharArray());
		HIVAAMutation pr67NMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 67, 'N');
		HIVAAMutation pr67XMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 67, 'X');
		HIVAAMutation pr67ADMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 67, "AD".toCharArray());
		HIVAAMutation pr67ADXMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 67, "AXD".toCharArray());
		HIVAAMutation eDiffN = pr67ADXMut;
		HIVAAMutation eDiffX = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 67, "ADN".toCharArray());
		HIVAAMutation eDiffAD = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 67, "XN".toCharArray());
		HIVAAMutation eDiffADX = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 67, 'N');
		assertEquals(eDiffN, pr67ANXDMut.subtractsBy(pr67NMut));
		assertEquals(eDiffX, pr67ANXDMut.subtractsBy(pr67XMut));
		assertEquals(eDiffAD, pr67ANXDMut.subtractsBy(pr67ADMut));
		assertEquals(eDiffADX, pr67ANXDMut.subtractsBy(pr67ADXMut));
	}

	@Test
	public void testSubtractsByEdgeCases() {
		HIVAAMutation pr68AMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, 'A');
		assertEquals(null, new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, 'A').subtractsBy(pr68AMut));
	}

	@Test(expected=IllegalArgumentException.class)
	public void testSubtractsByNull() {
		new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 67, 'A').subtractsBy((HIVAAMutation) null);
	}

	@Test(expected=IllegalArgumentException.class)
	public void testSubtractsByNotSamePos() {
		HIVAAMutation pr67AMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 67, 'A');
		new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, 'A').subtractsBy(pr67AMut);
	}

	@Test(expected=IllegalArgumentException.class)
	public void testSubtractsByNotSameGene() {
		HIVAAMutation rt67AMut = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 68, 'A');
		new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 67, 'A').subtractsBy(rt67AMut);
	}

	@Test(expected=IllegalArgumentException.class)
	public void testIntersectsWithNotSameGene() {
		new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, "AC".toCharArray())
			.intersectsWith(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 68, 'C'));
	}

	@Test(expected=IllegalArgumentException.class)
	public void testIntersectsWithNotSamePos() {
		new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, "AC".toCharArray())
			.intersectsWith(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 68, 'C'));
	}

	@Test
	public void testIntersectsWith() {
		assertEquals(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, 'C'),
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, "AC".toCharArray())
				.intersectsWith(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, "CD".toCharArray())));
	}

	@Test
	public void testSplit() {
		Set<HIVAAMutation> expected = new HashSet<>();
		expected.add(new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 151, 'L'));
		expected.add(new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 151, 'M'));
		expected.add(new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 151, 'Q'));
		assertEquals(expected, new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 151, "LMQ".toCharArray()).split());
		expected = new HashSet<>();
		expected.add(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, '-'));
		assertEquals(expected, new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, '-').split());
		expected = new HashSet<>();
		expected.add(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, '_'));
		assertEquals(expected, new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, '_').split());
		expected = new HashSet<>();
		expected.add(new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 1, 'L'));
		assertEquals(expected, new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 1, "PL".toCharArray()).split());
	}

	@Test
	public void testIsAtDrugResistancePosition() {
		HIVAAMutation mutIN151Major = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 151, 'L');
		assertTrue(mutIN151Major.isAtDrugResistancePosition());
		assertTrue(mutIN151Major.isAtDrugResistancePosition());

		HIVAAMutation mutIN150Other = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 150, 'L');
		assertFalse(mutIN150Other.isAtDrugResistancePosition());
		assertFalse(mutIN150Other.isAtDrugResistancePosition());

		HIVAAMutation mutIN151Accessory = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 151, 'A');
		assertTrue(mutIN151Accessory.isAtDrugResistancePosition());

		HIVAAMutation mutIN151Other = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 151, 'W');
		assertTrue(mutIN151Other.isAtDrugResistancePosition());

		HIVAAMutation mutIN151InsertionOther = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 151, "_L".toCharArray());
		assertTrue(mutIN151InsertionOther.isAtDrugResistancePosition());

		HIVAAMutation mutIN151DeletionOther = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 151, '-');
		assertTrue(mutIN151DeletionOther.isAtDrugResistancePosition());

		HIVAAMutation mutIN151MixtureWithMajor = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 151, "ALW".toCharArray());
		assertTrue(mutIN151MixtureWithMajor.isAtDrugResistancePosition());

		HIVAAMutation mutIN151MixtureWithoutMajor = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 151, "AW".toCharArray());
		assertTrue(mutIN151MixtureWithoutMajor.isAtDrugResistancePosition());

		HIVAAMutation mutRT75NRTI = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 75, 'I');
		assertTrue(mutRT75NRTI.isAtDrugResistancePosition());

		HIVAAMutation mutRT98NNRTI = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 98, 'G');
		assertTrue(mutRT98NNRTI.isAtDrugResistancePosition());

		HIVAAMutation mutRT99Other = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 99, 'G');
		assertFalse(mutRT99Other.isAtDrugResistancePosition());
	}

	@Test
	public void testGetDisplayAAs() {
		assertEquals("N", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, 'N').getDisplayAAs());
		assertEquals("X", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "ACDEFG".toCharArray(), 4).getDisplayAAs());
		assertEquals("ACDEFG", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "ACDEFG".toCharArray(), 7).getDisplayAAs());
	}

	@Test
	public void testGetAAs() {
		assertEquals("N", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, 'N').getAAs());
		assertEquals("ACDEFG", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "ACDEFG".toCharArray(), 4).getAAs());
		assertEquals("ACDEFG", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "ACDEFG".toCharArray(), 7).getAAs());
	}

	@Test
	public void testGetAAsWithoutConsensus() {
		assertEquals("N", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "KN".toCharArray()).getAAsWithoutReference());
	}

	@Test
	public void testGetTriplet() {
		// not support for an HIVMutation
		assertEquals("", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, 'A').getTriplet());
	}

	@Test
	public void testGetInsertedNAs() {
		// not support for an HIVMutation
		assertEquals("", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, '_').getInsertedNAs());
	}

	@Test
	public void testIsApobecDRM() {
		assertTrue(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 184, 'I').isApobecDRM());
		assertFalse(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 184, 'A').isApobecDRM());
		assertFalse(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 183, 'I').isApobecDRM());
	}

	@Test
	public void testIsApobecMut() {
		assertTrue(new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 25, "KN".toCharArray()).isApobecMutation());
		assertTrue(new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 25, "DN".toCharArray()).isApobecMutation());
		assertTrue(new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 25, "DKN".toCharArray()).isApobecMutation());
		assertFalse(new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 25, 'D').isApobecMutation());
		assertFalse(new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 25, 'A').isApobecMutation());
	}

	@Test
	public void testContainsSharedAA() {
		// contains shared consensus
		assertFalse(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "KN".toCharArray())
			.containsSharedAA(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "KA".toCharArray())));
		assertTrue(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "BKN".toCharArray())
			.containsSharedAA(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "NA".toCharArray())));
		assertFalse(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "BK".toCharArray())
			.containsSharedAA(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "NA".toCharArray())));
		assertFalse(
			new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 65, "NA".toCharArray())
			.containsSharedAA(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "NA".toCharArray())));
		assertFalse(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, "NA".toCharArray())
			.containsSharedAA(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "NA".toCharArray())));
		assertFalse(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, '*')
			.containsSharedAA(new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, '*')));
		assertFalse(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, '*')
			.containsSharedAA(Sets.newHashSet('*'), true));
		assertTrue(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, '*')
			.containsSharedAA(Sets.newHashSet('*'), false));
		assertFalse(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, "TD".toCharArray())
			.containsSharedAA(Sets.newHashSet('T'), true));
		assertTrue(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, "TD".toCharArray())
			.containsSharedAA(Sets.newHashSet('T'), false));
	}

	@Test
	public void testGetShortText() {
		assertEquals(
			"T69i", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, '_').getShortText());
		assertEquals(
			"D67d", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, '-').getShortText());
		assertEquals(
			"K65KA", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "KA".toCharArray()).getShortText());
		assertEquals(
			"K65KA", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "AK".toCharArray()).getShortText());
		assertEquals(
			"K65N", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, 'N').getShortText());
	}

	@Test
	public void testGetASIFormat() {
		assertEquals(
			"T69i", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, '_').getASIFormat());
		assertEquals(
			"D67d", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, '-').getASIFormat());
		assertEquals(
			"K65Z", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, '*').getASIFormat());
	}

	@Test
	public void testGetHIVDBFormat() {
		assertEquals(
			"69#", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, '_').getHIVDBFormat());
		assertEquals(
			"67~", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, '-').getHIVDBFormat());
		assertEquals(
			"65N", new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, 'N').getHIVDBFormat());
	}

	@Test
	public void testGetTypes() {
		HIVAAMutation mut1 = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 50, "VEF".toCharArray());
		HIVAAMutation mut2 = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 51, "VEF".toCharArray());
		List<MutationType> eTyoes1 = new ArrayList<>();
		eTyoes1.add(MutationType.Major);
		eTyoes1.add(MutationType.Accessory);
		List<MutationType>eTyoes2 = new ArrayList<>();
		eTyoes2.add(MutationType.Other);
		assertEquals(eTyoes1, mut1.getTypes());
		assertEquals(eTyoes1, mut1.getTypes()); // post instantiation of types
		assertEquals(eTyoes2, mut2.getTypes());
	}

	@Test
	public void testGetPrimaryType() {
		final HIVAAMutation majorMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 50, "VEF".toCharArray());
		final HIVAAMutation otherMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 51, "VEF".toCharArray());
		assertEquals(MutationType.Major, majorMut.getPrimaryType());
		assertEquals(MutationType.Other, otherMut.getPrimaryType());
	}

	@Test
	public void testEqualsAndHashCode() {
		final HIVAAMutation mut1 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, '_');
		final HIVAAMutation mut2 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, 'i');
		assertEquals(mut1, mut2);
		assertEquals(mut1, mut1);
		assertNotEquals(mut1, null);
		assertNotEquals(mut1, "T69_");
	}

	@Test
	public void testGetHumanFormat() {
		HIVAAMutation mut1 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "KN".toCharArray());
		HIVAAMutation mut2 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "NK".toCharArray());
		HIVAAMutation mut3 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 118, Sets.newHashSet('_'));
		HIVAAMutation mut4 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 118, Sets.newHashSet('#'));
		HIVAAMutation mut6 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, '-');
		HIVAAMutation mut8 = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 155, "S".toCharArray());
		HIVAAMutation mut9 = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 155, "NS".toCharArray());
		HIVAAMutation mut10 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 10, 'S');
		HIVAAMutation mut11 = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 10, 'S');
		HIVAAMutation mut12 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 215, "FIST".toCharArray());
		HIVAAMutation mut13 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 215, Sets.newHashSet('T', 'S', 'N', 'Y'));
		HIVAAMutation mut14 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 188, "YL*".toCharArray());
		HIVAAMutation mut15 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 188, '*');
		HIVAAMutation mut16 = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 263, "RKGY".toCharArray());
		HIVAAMutation mut17 = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 263, 'X');

		assertEquals("K65KN", mut1.getHumanFormat());
		assertEquals("K65KN", mut2.getHumanFormat());
		assertEquals(mut1, mut2);
		assertEquals(mut1.hashCode(), mut2.hashCode());

		assertEquals("V118Insertion", mut3.getHumanFormat());
		assertEquals("V118Insertion", mut4.getHumanFormat());
		assertNotEquals(mut1, mut3);
		assertNotEquals(mut2, mut4);
		assertEquals(mut3, mut4);
		// we consider insertions are the same
		assertEquals(mut3.hashCode(), mut4.hashCode());

		assertEquals("T69Deletion", mut6.getHumanFormat());

		assertEquals("N155S", mut8.getHumanFormat());
		assertEquals("N155NS", mut9.getHumanFormat());
		assertNotEquals(mut8, mut9);

		assertEquals("V10S", mut10.getHumanFormat());
		assertEquals("L10S", mut11.getHumanFormat());
		assertNotEquals(mut10, mut11);
		assertNotEquals(mut10.hashCode(), mut11.hashCode());

		assertEquals("T215TFIS", mut12.getHumanFormat());
		assertEquals("T215TNSY", mut13.getHumanFormat());
		assertNotEquals(mut12, mut13);
		assertNotEquals(mut12.hashCode(), mut13.hashCode());

		assertEquals("Y188Y*L", mut14.getHumanFormat());
		assertEquals("Y188*", mut15.getHumanFormat());
		assertNotEquals(mut14, mut15);
		assertNotEquals(mut14.hashCode(), mut15.hashCode());

		// TODO: shouldn't these two the same?
		assertEquals("R263RGKY", mut16.getHumanFormat());
		assertEquals("R263X", mut17.getHumanFormat());
		assertNotEquals(mut16, mut17);
		assertNotEquals(mut16.hashCode(), mut17.hashCode());
	}

	@Test
	public void testGetHumanFormatWithGene() {
		HIVAAMutation mut1 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "KN".toCharArray());
		HIVAAMutation mut2 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 65, "NK".toCharArray());
		HIVAAMutation mut3 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 118, "_".toCharArray());
		HIVAAMutation mut4 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 118, "#".toCharArray());
		HIVAAMutation mut6 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, "-".toCharArray());
		HIVAAMutation mut8 = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 155, "S".toCharArray());
		HIVAAMutation mut9 = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 155, "NS".toCharArray());
		HIVAAMutation mut10 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 10, "S".toCharArray());
		HIVAAMutation mut11 = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 10, "S".toCharArray());
		HIVAAMutation mut12 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 215, "FIST".toCharArray());
		HIVAAMutation mut13 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 215, "TSNY".toCharArray());
		HIVAAMutation mut14 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 188, "YL*".toCharArray());
		HIVAAMutation mut15 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 188, "*".toCharArray());
		HIVAAMutation mut16 = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 263, "RKGY".toCharArray());
		HIVAAMutation mut17 = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 263, "X".toCharArray());
		assertEquals("RT_K65KN", mut1.getHumanFormatWithGene());
		assertEquals("RT_K65KN", mut2.getHumanFormatWithGene());
		assertEquals("RT_V118Insertion", mut3.getHumanFormatWithGene());
		assertEquals("RT_V118Insertion", mut4.getHumanFormatWithGene());
		assertEquals("RT_T69Deletion", mut6.getHumanFormatWithGene());
		assertEquals("IN_N155S", mut8.getHumanFormatWithGene());
		assertEquals("IN_N155NS", mut9.getHumanFormatWithGene());
		assertEquals("RT_V10S", mut10.getHumanFormatWithGene());
		assertEquals("PR_L10S", mut11.getHumanFormatWithGene());
		assertEquals("RT_T215TFIS", mut12.getHumanFormatWithGene());
		assertEquals("RT_T215TNSY", mut13.getHumanFormatWithGene());
		assertEquals("RT_Y188Y*L", mut14.getHumanFormatWithGene());
		assertEquals("RT_Y188*", mut15.getHumanFormatWithGene());
		assertEquals("IN_R263RGKY", mut16.getHumanFormatWithGene());
		assertEquals("IN_R263X", mut17.getHumanFormatWithGene());
	}

	@Test
	public void testGetHumanFormatWithoutCons() {
		assertEquals(
			"69Insertion",
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, '#').getHumanFormatWithoutLeadingRef());
		assertEquals(
			"67Deletion",
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, '-').getHumanFormatWithoutLeadingRef());
	}

	@Test
	public void testIsIndel() {
		assertTrue(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, '_').isIndel());
		assertTrue(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, '-').isIndel());
		assertFalse(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, 'V').isIndel());
		assertTrue(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, new char[] {'i', 'V'}).isIndel());
	}

	@Test
	public void testIsMixture() {
		assertTrue(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 66, "BE".toCharArray()).isMixture());
		assertTrue(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 66, 'X').isMixture());
		assertFalse(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 66, 'V').isMixture());
		assertTrue(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, new char[] {'i', 'V'}).isMixture());
	}

	@Test
	public void testHasReference() {
		assertTrue(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, 'T').hasReference());
		assertTrue(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, "TV".toCharArray()).hasReference());
		assertFalse(
			new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, 'V').hasReference());
	}

	@Test
	public void testIsUnsequenced() {
		final HIVAAMutation mut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 1, 'X');
		// always false in HIVMutation
		assertFalse(mut.isUnsequenced());
	}

	@Test
	public void testGenePosition() {
		final HIVAAMutation mutPR68 = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, 'N');
		final HIVAAMutation mutRT67 = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, 'N');
		final HIVAAMutation mutIN155 = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 155, 'N');
		assertEquals(mutPR68.getGenePosition(), new HIVGenePosition(HIVGene.valueOf("HIV1PR"), 68));
		assertEquals(mutRT67.getGenePosition(), new HIVGenePosition(HIVGene.valueOf("HIV1RT"), 67));
		assertEquals(mutIN155.getGenePosition(), new HIVGenePosition(HIVGene.valueOf("HIV1IN"), 155));
	}

	@Test
	public void testIsInsertion() {
		final HIVAAMutation ins = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, '_');
		final HIVAAMutation insMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, "N_".toCharArray());
		final HIVAAMutation del = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, '-');
		final HIVAAMutation delIns = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, "_-".toCharArray());
		final HIVAAMutation mut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, 'N');
		assertTrue(ins.isInsertion());
		assertFalse(del.isInsertion());
		assertFalse(mut.isInsertion());
		assertTrue(insMut.isInsertion());
		assertTrue(delIns.isInsertion());
	}

	@Test
	public void testIsDeletion() {
		final HIVAAMutation del = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, '-');
		final HIVAAMutation delMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, "N-".toCharArray());
		final HIVAAMutation ins = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, '_');
		final HIVAAMutation delIns = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, "_-".toCharArray());
		final HIVAAMutation mut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, 'N');
		assertTrue(del.isDeletion());
		assertFalse(ins.isDeletion());
		assertFalse(mut.isDeletion());
		assertTrue(delMut.isDeletion());
		assertTrue(delIns.isDeletion());
	}

	@Test
	public void testHasStop() {
		final HIVAAMutation mut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 68, 'N');
		final HIVAAMutation stop = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 67, '*');
		final HIVAAMutation stopMut = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 155, "N*".toCharArray());
		assertFalse(mut.hasStop());
		assertTrue(stop.hasStop());
		assertTrue(stopMut.hasStop());
	}

	@Test
	public void testIsUnusual() {
		final HIVAAMutation unusualMut = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 1, 'A');
		final HIVAAMutation usualMut = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 1, 'L');
		final HIVAAMutation usualMuts = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 1, "LPS".toCharArray());
		final HIVAAMutation unusualMuts = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 1, "ACDEFG".toCharArray());
		final HIVAAMutation mixedMuts = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 75, "AILMVSTY".toCharArray());
		final HIVAAMutation unusualMutX = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 1, 'X');
		assertFalse(String.format("RT:%s should be usual", usualMut.toString()), usualMut.isUnusual());
		assertFalse(String.format("RT:%s should be usual", usualMuts.toString()), usualMuts.isUnusual());
		assertTrue(String.format("RT:%s should be unusual", unusualMut.toString()), unusualMut.isUnusual());
		assertTrue(String.format("RT:%s should be unusual", unusualMuts.toString()), unusualMuts.isUnusual());
		assertTrue(String.format("PR:%s should contain unusual mutation", mixedMuts.toString()), mixedMuts.isUnusual());
		assertTrue(String.format("RT:%s should be unusual", unusualMutX.toString()), unusualMutX.isUnusual());
	}

	@Test
	public void testIsSDRM() {
		final HIVAAMutation mut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 24, 'N');
		final HIVAAMutation sdrmMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 24, 'I');
		final HIVAAMutation mixedMuts = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 230, "LI".toCharArray());
		assertFalse(mut.isSDRM());
		assertTrue(sdrmMut.isSDRM());
		assertTrue(mixedMuts.isSDRM());
	}

	@Test
	public void testIsDRM() {
		final HIVAAMutation RT69ins = new HIVAAMutation(HIVGene.valueOf("HIV1RT"), 69, 'i');
		final HIVAAMutation IN263NK = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 263, "NK".toCharArray());
		assertTrue(RT69ins.isDRM());
		assertTrue(IN263NK.isDRM());
	}

	@Test
	public void testHasBDHVN() {
		final HIVAAMutation mut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 24, 'N');
		// always false
		assertFalse(mut.hasBDHVN());
	}

	@Test
	public void testIsAmbiguous() {
		final HIVAAMutation mut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 24, 'N');
		final HIVAAMutation xMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 24, 'X');
		final HIVAAMutation haMut = new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 24, "ABDEFGH".toCharArray());
		assertFalse(mut.isAmbiguous());
		assertTrue(xMut.isAmbiguous());
		assertTrue(haMut.isAmbiguous());
	}

	@Test
	public void testGetHighestMutPrevalance() {
		// Since we update prevalence data periodically, we
		// expects the following assertions to ultimately fail.
		// Hence we must manually update these assertions every time
		// we upload new prevalence data.
		final HIVAAMutation prevMut = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 45, 'G');
		final HIVAAMutation prevMuts = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 45, "HKQ".toCharArray());
		final HIVAAMutation prevMutZero = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 45, 'C');
		final HIVAAMutation prevMutsZero = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 45, "CDEFH".toCharArray());
		final HIVAAMutation prevMutsWCons = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 45, "LHKQ".toCharArray());
		final HIVAAMutation prevMutsWConsAndStop = new HIVAAMutation(HIVGene.valueOf("HIV1IN"), 45, "*LHKQ".toCharArray());
		assertEquals(0.03118, prevMut.getHighestMutPrevalence(), 1e-5);
		assertEquals(3.92463, prevMuts.getHighestMutPrevalence(), 1e-5);
		assertEquals(0.0, prevMutZero.getHighestMutPrevalence(), 1e-5);
		assertEquals(0.00445, prevMutsZero.getHighestMutPrevalence(), 1e-5);
		assertEquals(3.92463, prevMutsWCons.getHighestMutPrevalence(), 1e-5);
		assertEquals(3.92463, prevMutsWConsAndStop.getHighestMutPrevalence(), 1e-5);
		assertEquals(0.0, new HIVAAMutation(HIVGene.valueOf("HIV1PR"), 1, "PX".toCharArray()).getHighestMutPrevalence(), 1e-5);
	}
}
