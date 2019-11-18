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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.builder.HashCodeBuilder;

/**
 *
 * enum for Gene. Will replace the previous class-wrapper-around enum implementation now called Gene
 *
 */
public class HIVGene extends Gene<HIVStrain, HIVGene, HIVGenePosition> {

	private static final Map<String, HIVGene> singletons;
	private static final Character WILDCARD = '.';

	static {
		Map<String, HIVGene> _singletons = new HashMap<>();

		// LANL Consensus B
		_singletons.put(
			"HIV1PR",
			new HIVGene(
				HIVStrain.HIV1, HIVGeneEnum.PR, // 99 AAs
				"PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGI" +
				"GGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF",
				new Integer[] {},
				/* firstNA */2253));
		_singletons.put(
			"HIV1RT",
			new HIVGene(
				HIVStrain.HIV1, HIVGeneEnum.RT, // 560 AAs
				"PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKI" +
				"GPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGL" +
				"KKKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLP" +
				"QGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRT" +
				"KIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKD" +
				"SWTVNDIQKLVGKLNWASQIYAGIKVKQLCKLLRGTKALTEVIPLTEEAE" +
				"LELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLK" +
				"TGKYARMRGAHTNDVKQLTEAVQKIATESIVIWGKTPKFKLPIQKETWEA" +
				"WWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRET" +
				"KLGKAGYVTDRGRQKVVSLTDTTNQKTELQAIHLALQDSGLEVNIVTDSQ" +
				"YALGIIQAQPDKSESELVSQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDK" +
				"LVSAGIRKVL",
				new Integer[] {},
				/* firstNA */2550));
		_singletons.put(
			"HIV1IN",
			new HIVGene(
				HIVStrain.HIV1, HIVGeneEnum.IN, // 288 AAs
				"FLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAM" +
				"HGQVDCSPGIWQLDCTHLEGKIILVAVHVASGYIEAEVIPAETGQETAYF" +
				"LLKLAGRWPVKTIHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGV" +
				"VESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERI" +
				"VDIIATDIQTKELQKQITKIQNFRVYYRDSRDPLWKGPAKLLWKGEGAVV" +
				"IQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED",
				new Integer[] {},
				/* firstNA */4230));

		// ROD
		_singletons.put(
			"HIV2APR",
			new HIVGene(
				HIVStrain.HIV2A, HIVGeneEnum.PR, // 99 AAs
				"PQFSLWKRPVVTAYIEGQPVEVLLDTGADDSIVAGIELGNNYSPKIVGGI" +
				"GGFINTKEYKNVEIEVLNKKVRATIMTGDTPINIFGRNILTALGMSLNL",
				new Integer[] {},
				/* firstNA */2253)); // TODO: update firstNA to ROD
		_singletons.put(
			"HIV2ART",
			new HIVGene(
				HIVStrain.HIV2A, HIVGeneEnum.RT, // 559 AAs
				"PVAKVEPIKIMLKPGKDGPKLRQWPLTKEKIEALKEICEKMEKEGQLEEA" +
				"PPTNPYNTPTFAIKKKDKNKWRMLIDFRELNKVTQDFTEIQLGIPHPAGL" +
				"AKKRRITVLDVGDAYFSIPLHEDFRPYTAFTLPSVNNAEPGKRYIYKVLP" +
				"QGWKGSPAIFQHTMRQVLEPFRKANKDVIIIQYMDDILIASDRTDLEHDR" +
				"VVLQLKELLNGLGFSTPDEKFQKDPPYHWMGYELWPTKWKLQKIQLPQKE" +
				"IWTVNDIQKLVGVLNWAAQLYPGIKTKHLCRLIRGKMTLTEEVQWTELAE" +
				"AELEENRIILSQEQEGHYYQEEKELEATVQKDQENQWTYKIHQEEKILKV" +
				"GKYAKVKNTHTNGIRLLAQVVQKIGKEALVIWGRIPKFHLPVEREIWEQW" +
				"WDNYWQVTWIPDWDFVSTPPLVRLAFNLVGDPIPGAETFYTDGSCNRQSK" +
				"EGKAGYVTDRGKDKVKKLEQTTNQQAELEAFAMALTDSGPKVNIIVDSQY" +
				"VMGISASQPTESESKIVNQIIEEMIKKEAIYVAWVPAHKGIGGNQEVDHL" +
				"VSQGIRQVL",
				new Integer[] {
					/* insert one codon after 345 (RT346del) */ +345, 1
				},
				/* firstNA */2550)); // TODO: update firstNA to ROD
		_singletons.put(
			"HIV2AIN",
			new HIVGene(
				HIVStrain.HIV2A, HIVGeneEnum.IN, // 293 AAs
				"FLEKIEPAQEEHEKYHSNVKELSHKFGIPNLVARQIVNSCAQCQQKGEAI" +
				"HGQVNAELGTWQMDCTHLEGKIIIVAVHVASGFIEAEVIPQESGRQTALF" +
				"LLKLASRWPITHLHTDNGANFTSQEVKMVAWWIGIEQSFGVPYNPQSQGV" +
				"VEAMNHHLKNQISRIREQANTIETIVLMAIHCMNFKRRGGIGDMTPSERL" +
				"INMITTEQEIQFLQAKNSKLKDFRVYFREGRDQLWKGPGELLWKGEGAVL" +
				"VKVGTDIKIIPRRKAKIIRDYGGRQEMDSGSHLEGAREDGEMA",
				new Integer[] {
					/* delete two codons after 272 (IN272ins) */ -272, 2,
					/* delete three codons after 283 (IN283ins) */-283, 1,
					/* delete all codons after 288 */ -288, 0 },
				/* firstNA */4230)); // TODO: update firstNA to ROD

		// EHO
		_singletons.put(
			"HIV2BPR",
			new HIVGene(
				HIVStrain.HIV2B, HIVGeneEnum.PR, // 99 AAs
				"PQFSLWRRPVVKATIEGQSVEVLLDTGADDSIVAGIELGSNYTPKIVGGI" +
				"GGFINTNEYKNVEIEVVGKRVRATVMTGDTPINIFGRNILNSLGMTLNF",
				new Integer[] {},
				/* firstNA */2253)); // TODO: update firstNA to EHO
		_singletons.put(
			"HIV2BRT",
			new HIVGene(
				HIVStrain.HIV2B, HIVGeneEnum.RT, // 559 AAs
				"PVARIEPVKVQLKPEKDGPKIRQWPLSKEKILALKEICEKMEKEGQLEEA" +
				"PPTNPYNSPTFAIKKKDKNKWRMLIDFRELNKVTQEFTEVQLGIPHPAGL" +
				"ASKKRITVLDVGDAYFSVPLDPDFRQYTAFTLPAVNNAEPGKRYLYKVLP" +
				"QGWKGSPAIFQYTMAKVLDPFRKANNDVTIIQYMDDILVASDRSDLEHDR" +
				"VVSQLKELLNNMGFSTPEEKFQKDPPFKWMGYELWPKKWKLQKIQLPEKE" +
				"VWTVNDIQKLVGVLNWAAQLFPGIKTRHICKLIRGKMTLTEEVQWTELAE" +
				"AEFQENKIILEQEQEGSYYKEGVPLEATVQKNLANQWTYKIHQGDKILKV" +
				"GKYAKVKNTHTNGVRLLAHVVQKIGKEALVIWGEIPMFHLPVERETWDQW" +
				"WTDYWQVTWIPEWDFVSTPPLIRLAYNLVKDPLEGVETYYTDGSCNKASK" +
				"EGKAGYVTDRGKDKVKPLEQTTNQQAELEAFALALQDSGPQVNIIVDSQY" +
				"VMGIVAAQPTETESPIVREIIEEMIKKEKIYVGWVPAHKGLGGNQEVDHL" +
				"VSQGIRQIL",
				new Integer[] {
					/* insert one codon after 345 (RT346del) */ +345, 1
				},
				/* firstNA */2550)); // TODO: update firstNA to EHO
		_singletons.put(
			"HIV2BIN",
			new HIVGene(
				HIVStrain.HIV2B, HIVGeneEnum.IN, // 296 AAs
				"FLEKIEPAQEEHEKYHNNVKELVHKFGIPQLVARQIVNSCDKCQQKGEAI" +
				"HGQVNSELGTWQMDCTHLEGKVIIVAVHVASGFIEAEVIPQETGRQTALF" +
				"LLKLASRWPITHLHTDNGANFTSQDVKMAAWWIGIEQTFGVPYNPESQGV" +
				"VEAMNHHLKNQIDRIRDQAVSIETVVLMATHCMNFKRRGGIGDMTPAERI" +
				"VNMITTEQEIQFLQTKNLKFQNFRVYYREGRDQLWKGPGDLLWKGEGAVI" +
				"IKVGTEIKVIPRRKAKIIRNYGGGKELDCSADVEDTMQAREVAQSN",
				new Integer[] {
					/* delete two codons after 272 (IN272ins) */ -272, 2,
					/* delete three codons after 283 (IN283ins) */-283, 1,
					/* delete all codons after 288 */ -288, 0 },
				/* firstNA */4230)); // TODO: update firstNA to EHO


		singletons = Collections.unmodifiableMap(_singletons);
	}

	private final HIVStrain strain;
	private final HIVGeneEnum geneEnum;
	private final String reference;
	private final Integer[] alignmentAdjustment;
	private final int firstNA;

	protected HIVGene(
		HIVStrain strain, HIVGeneEnum gene, String reference,
		Integer[] alignmentAdjustment, int firstNA
	) {
		this.strain = strain;
		this.geneEnum = gene;
		this.reference = reference;
		this.alignmentAdjustment = alignmentAdjustment;
		this.firstNA = firstNA;
	}
	public static HIVGene valueOf(String strainText, String geneText) {
		return singletons.get(String.format("%s%s", strainText, geneText));
	}

	public static HIVGene valueOf(HIVStrain strain, String geneText) {
		return singletons.get(String.format("%s%s", strain, geneText));
	}

	public static HIVGene valueOf(HIVStrain strain, HIVGeneEnum gene) {
		return singletons.get(String.format("%s%s", strain, gene));
	}

	public static HIVGene valueOf(String strainGeneText) {
		return singletons.get(strainGeneText);
	}

	public static HIVGene[] values(String strainText) {
		return new HIVGene[] {
			singletons.get(strainText + "PR"),
			singletons.get(strainText + "RT"),
			singletons.get(strainText + "IN")
		};
	}

	public static HIVGene[] values(HIVStrain strain) {
		return new HIVGene[] {
			singletons.get(strain + "PR"),
			singletons.get(strain + "RT"),
			singletons.get(strain + "IN")
		};
	}

	/**
	 * Get the drug classes associated with a gene
	 */
	@Override
	public List<HIVDrugClass> getDrugClasses() {
		List<HIVDrugClass> drugClasses = new ArrayList<>();
		switch(geneEnum) {
		case RT:
			drugClasses.add(HIVDrugClass.NRTI);
			if (strain == HIVStrain.HIV1) {
				drugClasses.add(HIVDrugClass.NNRTI);
			}
			break;
		case PR:
			drugClasses.add(HIVDrugClass.PI);
			break;
		default: // case IN:
			drugClasses.add(HIVDrugClass.INSTI);
			break;
		}
		return drugClasses;
	}

	@Override
	public HIVGenePosition getGenePosition(Integer pos) {
		return new HIVGenePosition(this, pos);
	}

	/**
	 * Get the mutation types associated with a gene
	 */
	@Override
	public List<MutType> getMutationTypes() {
		List<MutType> mutTypes = new ArrayList<>();
		switch(geneEnum) {
		case RT:
			if (strain == HIVStrain.HIV1) {
				mutTypes.add(MutType.NRTI);
				mutTypes.add(MutType.NNRTI);
			} else {
				mutTypes.add(MutType.Major);
				mutTypes.add(MutType.Accessory);
			}
			mutTypes.add(MutType.Other);
			break;
		case PR:
			mutTypes.add(MutType.Major);
			mutTypes.add(MutType.Accessory);
			mutTypes.add(MutType.Other);
			break;
		default: // case IN:
			mutTypes.add(MutType.Major);
			mutTypes.add(MutType.Accessory);
			mutTypes.add(MutType.Other);
			break;
		}
		return mutTypes;
	}

	/**
	 * Get the scored mutation types associated with a gene
	 */
	@Override
	public List<MutType> getScoredMutTypes() {
		List<MutType> mutTypes = new ArrayList<>();
		switch(geneEnum) {
		case RT:
			if (strain == HIVStrain.HIV1) {
				mutTypes.add(MutType.NRTI);
				mutTypes.add(MutType.NNRTI);
			} else {
				mutTypes.add(MutType.Major);
				mutTypes.add(MutType.Accessory);
			}
			break;
		case PR:
			mutTypes.add(MutType.Major);
			mutTypes.add(MutType.Accessory);
			break;
		default: // case IN:
			mutTypes.add(MutType.Major);
			mutTypes.add(MutType.Accessory);
			break;
		}
		return mutTypes;
	}

	@Override
	public int getLength() {
		return this.reference.length();
	}

	@Override
	public int getFirstNA() {
		return this.firstNA;
	}

	@Override
	public int getNASize() {
		return this.reference.length() * 3;
	}

	/**
	 * Get the reference amino acid (AA) at a position in a gene
	 * Indexed starting from 1
	 * @param pos
	 * @return the AA at the submitted position
	 */
	@Override
	public String getReference(int pos) {
		return getRefChar(pos).toString();
	}

	@Override
	public Character getRefChar(int pos) {
		return this.reference.charAt(pos - 1);
	}

	@Override
	public String getReference(int pos, int length) {
		return this.reference.substring(pos - 1, pos - 1 + length);
	}

	@Override
	public String getReference() {
		return this.reference;
	}

	@Override
	public HIVStrain getStrain() {
		return this.strain;
	}

	public HIVGeneEnum getGeneEnum() {
		return this.geneEnum;
	}

	@Override
	public String getNameWithStrain() {
		return String.format("%s%s", strain, geneEnum);
	}

	@Override
	public String getName() {
		return this.geneEnum.name();
	}

	/**
	 * Adjust AA Alignment with HXB2 reference
	 *
	 * For RT and IN of HIV-2 virus, there are deletions and insertions at
	 * non-DRM positions comparing to HXB reference. This adjustment allows
	 * to "re-fit" the HIV-2 sequences into HXB2 numbering system. So these
	 * sequences can be compatible with the HIV-1 dedicated data structure.
	 *
	 * @param aaseq
	 * @param firstAA
	 * @param lastAA
	 * @return String of adjusted AA alignment in full gene length
	 */
	public String adjustAAAlignment(String aaseq, int firstAA, int lastAA) {
		int numPrefixAAs = firstAA - 1;
		int numSuffixAAs = this.reference.length() - aaseq.length() - numPrefixAAs;
		aaseq =
			StringUtils.repeat(WILDCARD, numPrefixAAs) +
			aaseq +
			StringUtils.repeat(WILDCARD, numSuffixAAs);
		if (alignmentAdjustment.length == 0) {
			return aaseq;
		}
		for (int i = 0; i < alignmentAdjustment.length; i ++) {
			int pos = alignmentAdjustment[i];
			int aaSize = alignmentAdjustment[++ i];
			if (pos > 0) {
				if (aaSize <= 0) {
					throw new RuntimeException(String.format(
						"unable to add %s AA(s) to aaseq (AlignmentAdjustment)", aaSize
					));
				}
				// HXB2 deletion, add AA(s) to aaseq
				aaseq =
					aaseq.substring(0, pos) +
					StringUtils.repeat(WILDCARD, aaSize) +
					aaseq.substring(pos);
			}
			else {
				// HXB2 insertion, delete AA(s) from aaseq
				pos = -pos;
				if (aaSize > 0) {
					aaseq =
						aaseq.substring(0, pos) +
						aaseq.substring(pos + aaSize);
				}
				else if (aaSize == 0) {
					// delete anything after
					aaseq = aaseq.substring(0, pos);
				}
				else {
					throw new RuntimeException(String.format(
						"unable to remove %s AA(s) from aaseq (AlignmentAdjustment)", aaSize
					));
				}
			}
		}
		return aaseq;
	}

	/**
	 * Adjust NA Alignment with HXB2 reference
	 *
	 * For RT and IN of HIV-2 virus, there are deletions and insertions at
	 * non-DRM positions comparing to HXB reference. This adjustment allows
	 * to "re-fit" the HIV-2 sequences into HXB2 numbering system. So these
	 * sequences can be compatible with the HIV-1 dedicated data structure.
	 *
	 * @param naseq
	 * @param firstAA
	 * @param lastAA
	 * @return String of adjusted NA alignment in full gene length
	 */
	public String adjustNAAlignment(String naseq, int firstAA, int lastAA) {
		int numPrefixNAs = (firstAA - 1) * 3;
		int numSuffixNAs = getNASize() - naseq.length() - numPrefixNAs;
		naseq =
			StringUtils.repeat(WILDCARD, numPrefixNAs) +
			naseq +
			StringUtils.repeat(WILDCARD, numSuffixNAs);
		if (alignmentAdjustment.length == 0) {
			return naseq;
		}
		for (int i = 0; i < alignmentAdjustment.length; i ++) {
			int pos = alignmentAdjustment[i];
			int aaSize = alignmentAdjustment[++ i];
			if (pos > 0) {
				if (aaSize <= 0) {
					throw new RuntimeException(String.format(
						"unable to add %s codon(s) to naseq (AlignmentAdjustment)", aaSize
					));
				}
				// HXB2 deletion, add codon(s) to naseq
				naseq =
					naseq.substring(0, pos * 3) +
					StringUtils.repeat(WILDCARD, 3 * aaSize) +
					naseq.substring(pos * 3);
			}
			else {
				// HXB2 insertion, delete codon(s) from naseq
				pos = -pos;
				if (aaSize > 0) {
					naseq =
						naseq.substring(0, pos * 3) +
						naseq.substring((pos + aaSize) * 3);
				}
				else if (aaSize == 0) {
					// delete anything after
					naseq = naseq.substring(0, pos * 3);
				}
				else {
					throw new RuntimeException(String.format(
						"unable to remove %s codon(s) from naseq (AlignmentAdjustment)", aaSize
					));
				}
			}
		}
		return naseq;
	}

	@Override
	public int compareTo(HIVGene o) {
		if (o == null) throw new NullPointerException("Null is incomprable.");
		HIVGene oo = (HIVGene) o;
		int cmp = strain.compareTo(oo.strain) * 3;
		cmp += geneEnum.compareTo(oo.geneEnum);
		return cmp;
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(23564345, 437136943)
			.append(strain).append(geneEnum).toHashCode();
	}

	@Override
	public boolean equals(Object o) {
		// singleton only has one reference point
		return this == o;
	}

}
