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

import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.Test;

import com.google.common.base.Strings;

import edu.stanford.hivdb.mutations.MutationType;

public class GeneTest {

	@Test
	public void testValues() {
		assertArrayEquals(HIVGene.values("HIV1"), new HIVGene[] {
			HIVGene.valueOf("HIV1PR"),
			HIVGene.valueOf("HIV1RT"),
			HIVGene.valueOf("HIV1IN")
		});
		assertArrayEquals(HIVGene.values(HIVStrain.HIV2B), new HIVGene[] {
			HIVGene.valueOf("HIV2BPR"),
			HIVGene.valueOf("HIV2BRT"),
			HIVGene.valueOf("HIV2BIN")
		});
		assertArrayEquals(HIVGene.values("HIV-1"), new HIVGene[] {
			null,
			null,
			null
		});
	}

	@Test
	public void testValueOf() {
		assertEquals(
			HIVGene.valueOf("HIV1", "PR"), HIVGene.valueOf(HIVStrain.HIV1, HIVAbstractGene.PR));
		assertEquals(
			HIVGene.valueOf("HIV2BRT"), HIVGene.valueOf(HIVStrain.HIV2B, HIVAbstractGene.RT));
		assertNotEquals(
			HIVGene.valueOf("HIV2ART"), HIVGene.valueOf(HIVStrain.HIV2B, HIVAbstractGene.RT));
		assertEquals(HIVGene.valueOf("HIV2CRT"), null);
		assertEquals(HIVGene.valueOf("HIV2A", "Pol"), null);
		assertEquals(HIVGene.valueOf(HIVStrain.HIV1, "Pol"), null);
	}

	@Test
	public void testGetDrugClasses() {
		assertEquals(
			Arrays.asList(new HIVDrugClass[] {HIVDrugClass.PI}),
			HIVGene.valueOf("HIV1PR").getDrugClasses());
		assertEquals(
			Arrays.asList(new HIVDrugClass[] {HIVDrugClass.NRTI, HIVDrugClass.NNRTI}),
			HIVGene.valueOf("HIV1RT").getDrugClasses());
		assertEquals(
			Arrays.asList(new HIVDrugClass[] {HIVDrugClass.NRTI}),
			HIVGene.valueOf("HIV2ART").getDrugClasses());
		assertEquals(
			Arrays.asList(new HIVDrugClass[] {HIVDrugClass.NRTI}),
			HIVGene.valueOf("HIV2BRT").getDrugClasses());
		assertEquals(
			Arrays.asList(new HIVDrugClass[] {HIVDrugClass.INSTI}),
			HIVGene.valueOf("HIV1IN").getDrugClasses());
	}

	@Test
	public void testGetMutationTypes() {
		assertEquals(
			Arrays.asList(new MutationType[] {MutationType.Major, MutationType.Accessory, MutationType.Other}),
			HIVGene.valueOf("HIV1PR").getMutationTypes());
		assertEquals(
			Arrays.asList(new MutationType[] {MutationType.NRTI, MutationType.NNRTI, MutationType.Other}),
			HIVGene.valueOf("HIV1RT").getMutationTypes());
		assertEquals(
			Arrays.asList(new MutationType[] {MutationType.Major, MutationType.Accessory, MutationType.Other}),
			HIVGene.valueOf("HIV2ART").getMutationTypes());
		assertEquals(
			Arrays.asList(new MutationType[] {MutationType.Major, MutationType.Accessory, MutationType.Other}),
			HIVGene.valueOf("HIV2BRT").getMutationTypes());
		assertEquals(
			Arrays.asList(new MutationType[] {MutationType.Major, MutationType.Accessory, MutationType.Other}),
			HIVGene.valueOf("HIV1IN").getMutationTypes());
	}

	@Test
	public void testGetScoredMutTypes() {
		assertEquals(
			Arrays.asList(new MutationType[] {MutationType.Major, MutationType.Accessory}),
			HIVGene.valueOf("HIV1PR").getScoredMutTypes());
		assertEquals(
			Arrays.asList(new MutationType[] {MutationType.NRTI, MutationType.NNRTI}),
			HIVGene.valueOf("HIV1RT").getScoredMutTypes());
		assertEquals(
			Arrays.asList(new MutationType[] {MutationType.Major, MutationType.Accessory}),
			HIVGene.valueOf("HIV2ART").getScoredMutTypes());
		assertEquals(
			Arrays.asList(new MutationType[] {MutationType.Major, MutationType.Accessory}),
			HIVGene.valueOf("HIV2BRT").getScoredMutTypes());
		assertEquals(
			Arrays.asList(new MutationType[] {MutationType.Major, MutationType.Accessory}),
			HIVGene.valueOf("HIV1IN").getScoredMutTypes());
	}

	@Test
	public void testGetLength() {
		assertEquals(99, HIVGene.valueOf("HIV1PR").getLength());
		assertEquals(560, HIVGene.valueOf("HIV1RT").getLength());
		assertEquals(288, HIVGene.valueOf("HIV1IN").getLength());
	}

	@Test
	public void testGetConsensus() {
		assertEquals("T", HIVGene.valueOf("HIV1PR").getReference(4));
		assertEquals("IDK", HIVGene.valueOf("HIV1IN").getReference(5, 3));
		assertEquals(560, HIVGene.valueOf("HIV1RT").getReference().length());
	}

	@Test
	public void testGetFirstNA() {
		assertEquals(2253, HIVGene.valueOf("HIV1PR").getFirstNA());
		assertEquals(2550, HIVGene.valueOf("HIV1RT").getFirstNA());
		assertEquals(4230, HIVGene.valueOf("HIV1IN").getFirstNA());
	}

	@Test
	public void testGetNASize() {
		assertEquals(297, HIVGene.valueOf("HIV1PR").getNASize());
		assertEquals(1680, HIVGene.valueOf("HIV1RT").getNASize());
		assertEquals(864, HIVGene.valueOf("HIV1IN").getNASize());
		assertEquals(297, HIVGene.valueOf("HIV2APR").getNASize());
		assertEquals(1677, HIVGene.valueOf("HIV2ART").getNASize());
		assertEquals(879, HIVGene.valueOf("HIV2AIN").getNASize());
		assertEquals(297, HIVGene.valueOf("HIV2BPR").getNASize());
		assertEquals(1677, HIVGene.valueOf("HIV2BRT").getNASize());
		assertEquals(888, HIVGene.valueOf("HIV2BIN").getNASize());
	}

	@Test
	public void testGetStrain() {
		assertEquals(HIVStrain.HIV2A, HIVGene.valueOf("HIV2APR").getStrain());
		assertEquals(HIVStrain.HIV2B, HIVGene.valueOf("HIV2BRT").getStrain());
		assertEquals(HIVStrain.HIV1, HIVGene.valueOf("HIV1IN").getStrain());
	}

	@Test
	public void testGetGeneEnum() {
		assertEquals(HIVAbstractGene.PR, HIVGene.valueOf("HIV2APR").getAbstractGene());
		assertEquals(HIVAbstractGene.RT, HIVGene.valueOf("HIV2BRT").getAbstractGene());
		assertEquals(HIVAbstractGene.IN, HIVGene.valueOf("HIV1IN").getAbstractGene());
	}

	@Test
	public void testGetNameWithStrain() {
		assertEquals("HIV1PR", HIVGene.values("HIV1")[0].getNameWithStrain());
		assertEquals("HIV2ART", HIVGene.values("HIV2A")[1].getNameWithStrain());
		assertEquals("HIV2BIN", HIVGene.values("HIV2B")[2].getNameWithStrain());
	}

	@Test
	public void testGetName() {
		assertEquals("PR", HIVGene.values("HIV1")[0].getName());
		assertEquals("RT", HIVGene.values("HIV2A")[1].getName());
		assertEquals("IN", HIVGene.values("HIV2B")[2].getName());
	}

	@Test
	public void testAdjustAAAlignment() {
		assertEquals(
			Strings.repeat(".", 4) + "ABCDEF" + Strings.repeat(".", 89),
			HIVGene.valueOf("HIV1PR").adjustAAAlignment("ABCDEF", 5, 10));
		assertEquals(
			Strings.repeat(".", 342) + "FEDCBA" + Strings.repeat(".", 212),
			HIVGene.valueOf("HIV1RT").adjustAAAlignment("FEDCBA", 343, 348));
		assertEquals(
			// RT346 Deletion
			Strings.repeat(".", 342) + "FED.CBA" + Strings.repeat(".", 211),
			HIVGene.valueOf("HIV2ART").adjustAAAlignment("FEDCBA", 343, 348));
		assertEquals(
			// RT346 Deletion
			Strings.repeat(".", 342) + "FED.CBA" + Strings.repeat(".", 211),
			HIVGene.valueOf("HIV2BRT").adjustAAAlignment("FEDCBA", 343, 348));
		assertEquals(
			// IN272 Insertion
			Strings.repeat(".", 270) + "ABEF" + Strings.repeat(".", 14),
			HIVGene.valueOf("HIV2AIN").adjustAAAlignment("ABCDEF", 271, 276));
		assertEquals(
			// IN283 Insertion + IN272 two AAs shift
			Strings.repeat(".", 278) + "ABCDEG....",
			HIVGene.valueOf("HIV2AIN").adjustAAAlignment("ABCDEFG", 281, 287));
		assertEquals(
			// IN after 288 (IN272 + IN283 three AAs shift)
			Strings.repeat(".", 283) + "ABCDE",
			HIVGene.valueOf("HIV2AIN").adjustAAAlignment("ABCDEFG", 287, 293));
	}

	/*@Test(expected = RuntimeException.class)
	public void testAdjustAAAlignmentWithException1() {
		HIVGene fakePRGene = new HIVGene(
			HIVStrain.HIV2A, HIVAbstractGene.PR,
			"PQFSLWRRPVVKATIEGQSVEVLLDTGADDSIVAGIELGSNYTPKIVGGI" +
			"GGFINTNEYKNVEIEVVGKRVRATVMTGDTPINIFGRNILNSLGMTLNF",
			new Integer[] {
				+22, -1
			}, 2253);
		fakePRGene.adjustAAAlignment("ABCDEF", 20, 25);
	}

	@Test(expected = RuntimeException.class)
	public void testAdjustAAAlignmentWithException2() {
		HIVGene fakePRGene = new HIVGene(
			HIVStrain.HIV2A, HIVAbstractGene.PR,
			"PQFSLWRRPVVKATIEGQSVEVLLDTGADDSIVAGIELGSNYTPKIVGGI" +
			"GGFINTNEYKNVEIEVVGKRVRATVMTGDTPINIFGRNILNSLGMTLNF",
			new Integer[] {
				-22, -1
			}, 2253);
		fakePRGene.adjustAAAlignment("ABCDEF", 20, 25);
	}*/

	@Test
	public void testAdjustNAAlignment() {
		assertEquals(
			Strings.repeat("...", 4) + "AAABBBCCCDDDEEEFFF" + Strings.repeat("...", 89),
			HIVGene.valueOf("HIV1PR").adjustNAAlignment("AAABBBCCCDDDEEEFFF", 5, 10));
		assertEquals(
			Strings.repeat("...", 342) + "FFFEEEDDDCCCBBBAAA" + Strings.repeat("...", 212),
			HIVGene.valueOf("HIV1RT").adjustNAAlignment("FFFEEEDDDCCCBBBAAA", 343, 348));
		assertEquals(
			// RT346 Deletion
			Strings.repeat("...", 342) + "FFFEEEDDD...CCCBBBAAA" + Strings.repeat("...", 211),
			HIVGene.valueOf("HIV2ART").adjustNAAlignment("FFFEEEDDDCCCBBBAAA", 343, 348));
		assertEquals(
			// RT346 Deletion
			Strings.repeat("...", 342) + "FFFEEEDDD...CCCBBBAAA" + Strings.repeat("...", 211),
			HIVGene.valueOf("HIV2BRT").adjustNAAlignment("FFFEEEDDDCCCBBBAAA", 343, 348));
		assertEquals(
			// IN272 Insertion
			Strings.repeat("...", 270) + "AAABBBEEEFFF" + Strings.repeat("...", 14),
			HIVGene.valueOf("HIV2AIN").adjustNAAlignment("AAABBBCCCDDDEEEFFF", 271, 276));
		assertEquals(
			// IN283 Insertion + IN272 two AAs shift
			Strings.repeat("...", 278) + "AAABBBCCCDDDEEEGGG............",
			HIVGene.valueOf("HIV2AIN").adjustNAAlignment("AAABBBCCCDDDEEEFFFGGG", 281, 287));
		assertEquals(
			// IN after 288 (IN272 + IN283 three AAs shift)
			Strings.repeat("...", 283) + "AAABBBCCCDDDEEE",
			HIVGene.valueOf("HIV2AIN").adjustNAAlignment("AAABBBCCCDDDEEEFFFGGG", 287, 293));
	}

	/*@Test(expected = RuntimeException.class)
	public void testAdjustNAAlignmentWithException1() {
		HIVGene fakePRGene = new HIVGene(
			HIVStrain.HIV2A, HIVAbstractGene.PR,
			"PQFSLWRRPVVKATIEGQSVEVLLDTGADDSIVAGIELGSNYTPKIVGGI" +
			"GGFINTNEYKNVEIEVVGKRVRATVMTGDTPINIFGRNILNSLGMTLNF",
			new Integer[] {
				+22, -1
			}, 2253);
		fakePRGene.adjustNAAlignment("ABCDEF", 20, 25);
	}

	@Test(expected = RuntimeException.class)
	public void testAdjustNAAlignmentWithException2() {
		HIVGene fakePRGene = new HIVGene(
			HIVStrain.HIV2A, HIVAbstractGene.PR,
			"PQFSLWRRPVVKATIEGQSVEVLLDTGADDSIVAGIELGSNYTPKIVGGI" +
			"GGFINTNEYKNVEIEVVGKRVRATVMTGDTPINIFGRNILNSLGMTLNF",
			new Integer[] {
				-22, -1
			}, 2253);
		fakePRGene.adjustNAAlignment("ABCDEF", 20, 25);
	}*/

	@Test
	public void testCompareTo() {
		assertEquals(0, HIVGene.valueOf("HIV1", "PR").compareTo(HIVGene.valueOf("HIV1PR")));
		assertEquals(-1, HIVGene.valueOf("HIV1", "PR").compareTo(HIVGene.valueOf("HIV1RT")));
		assertEquals(2, HIVGene.valueOf("HIV2AIN").compareTo(HIVGene.valueOf("HIV2APR")));
		assertEquals(-1, HIVGene.valueOf("HIV1IN").compareTo(HIVGene.valueOf("HIV2APR")));
		assertEquals(-2, HIVGene.valueOf("HIV1IN").compareTo(HIVGene.valueOf("HIV2ART")));
		assertEquals(-3, HIVGene.valueOf("HIV1IN").compareTo(HIVGene.valueOf("HIV2AIN")));
		assertEquals(-4, HIVGene.valueOf("HIV1IN").compareTo(HIVGene.valueOf("HIV2BPR")));
		assertEquals(-5, HIVGene.valueOf("HIV1RT").compareTo(HIVGene.valueOf("HIV2BPR")));
		assertEquals(4, HIVGene.valueOf("HIV2ART").compareTo(HIVGene.valueOf("HIV1PR")));
	}

	@Test(expected = NullPointerException.class)
	public void testCompareToNull() {
		HIVGene.valueOf("HIV1RT").compareTo(null);
	}

	@Test
	public void testHashCode() {
		assertEquals(HIVGene.values("HIV1")[0].hashCode(), HIVGene.valueOf("HIV1PR").hashCode());
	}

}
