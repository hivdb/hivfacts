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

package edu.stanford.hivdb.hivfacts.extras;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.lang3.StringUtils;

import edu.stanford.hivdb.utilities.TSV;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.hivfacts.HIV;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.mutations.FrameShift;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.sequences.AlignedGeneSeq;
import edu.stanford.hivdb.sequences.AlignedSequence;
import edu.stanford.hivdb.utilities.NumberFormats;


public class TabularSequenceSummary {
	private static final String[] headerFields = {
			"Sequence Name", "Genes", "PR Start", "PR End", "RT Start",
			"RT End", "IN Start", "IN End", "Subtype (%)", "Pcnt Mix",
			"PR Major", "PR Accessory", "PR Other", "NRTI Major", "NRTI Accessory",
			"NNRTI Major", "NNRTI Accessory", "RT Other", "IN Major", 
			"IN Accessory", "IN Other", "PR SDRMs", "RT SDRMs", "IN SDRMs",
			"PI TSMs", "NRTI TSMs", "NNRTI TSMs", "INSTI TSMs",
			"Num Frame Shifts", "Frame Shifts",	"Num Insertions", "Insertions",
			"Num Deletions", "Deletions", "Num Stop Codons", "StopCodons",
			"Num BDHVN", "BDHVN", "Num Apobec Mutations", "Apobec Mutations",
			"Num Unusual Mutations", "UnusualMutations"};
	private List<List<String>> sequenceRows = new ArrayList<>();
	private Map<String, Map<String, String>> tabularResults = new HashMap<>();

	/**
	 *
	 */
	public TabularSequenceSummary (List<AlignedSequence<HIV>> overallResults) {

		for (AlignedSequence<HIV> alignedSeq : overallResults) {
			List<String> sequenceRecord = new ArrayList<>();

			List<Gene<HIV>> geneList = alignedSeq.getAvailableGenes();
			MutationSet<HIV> seqMutations = alignedSeq.getMutations();
			String seqName = alignedSeq.getInputSequence().getHeader();

			tabularResults.put(seqName, new HashMap<String, String>());
			String genes = StringUtils.join(
				geneList.stream().map(g -> g.getName()).toArray(), ",");

			// sequenceName
			sequenceRecord.add(seqName);

			// Genes
			sequenceRecord.add(genes);

			// PRStart, PREnd, RTStart, RTEnd, INStart, INEnd
			sequenceRecord.addAll(determineGeneBoundaries(alignedSeq));

			// Subtype(%)
			sequenceRecord.add(determineSubtype(alignedSeq));

			// PcntMix
			sequenceRecord.add(
				NumberFormats.prettyDecimalAsString(alignedSeq.getMixturePcnt()));

			// PRMajor, PRAccessory, PROther,
			// NRTIMajor, NRTIAccessory, NNRTIMajor, NNRTIAccessory, RTOther,
			// INMajor, INAccessory, INOther
			sequenceRecord.addAll(determineMutLists(alignedSeq));

		  	// PRSDRMs, RTSDRMs
			sequenceRecord.addAll(determineSdrms(alignedSeq));

		  	// PI-TSMs, NRTI-TSMs, NNRTI-TSMs, INSTI-TSMs
			sequenceRecord.addAll(determineNonDrmTsms(alignedSeq));

		  	// NumFS, FrameShifts
			sequenceRecord.addAll(determineFrameShiftText(alignedSeq));
			// NumIns, Insertions
			sequenceRecord.addAll(determineSeqInsertions(seqMutations));
			// NumDel, Deletions
			sequenceRecord.addAll(determineSeqDeletions(seqMutations));
			// NumStops, StopCodons
			sequenceRecord.addAll(determineSeqStopCodons(seqMutations));
			// NumBDHVN, BDHVN
			sequenceRecord.addAll(determineSeqBDHVN(seqMutations));
			// NumApobec, ApobecMuts
			sequenceRecord.addAll(determineApobecFields(seqMutations));
			// NumUnusual, UnusualMuts
			sequenceRecord.addAll(determineSeqUnusualMuts(seqMutations));
			sequenceRows.add(sequenceRecord);

			for (int i=0; i<headerFields.length; i++) {
				String field = headerFields[i];
				String dataItem = sequenceRecord.get(i);
				tabularResults.get(seqName).put(field, dataItem);
			}
		}
	}

	@Override
	public String toString() {
		return TSV.dumps(headerFields, sequenceRows);
	}

	public String getHeader() {
		return TSV.dumpsHeader(headerFields);
	}

	public String getBody() {
		return TSV.dumpsBody(sequenceRows);
	}

	public Map<String, Map<String, String>> getTable() { return tabularResults; }
	public String[] getHeaderFields() { return headerFields; }

	private static String determineSubtype(AlignedSequence<HIV> alignedSeq) {
		return alignedSeq.getGenotypeText();
	}

	private static List<String> mutationListToTabularResult(MutationSet<HIV> mutations) {
		List<String> fields = new ArrayList<>();
		String text = "None";
		if (mutations.size() > 0) {
			text = mutations.join();
		}
		fields.add("" + mutations.size());
		fields.add(text);
		return fields;
	}

	private static List<String> determineApobecFields(MutationSet<HIV> mutations) {
		MutationSet<HIV> apobecMuts = mutations.getApobecMutations();
		return mutationListToTabularResult(apobecMuts);
	}

	private static List<String> determineSeqUnusualMuts(MutationSet<HIV> mutations) {
		MutationSet<HIV> unusualMuts = mutations.getUnusualMutations();
		return mutationListToTabularResult(unusualMuts);
	}

	// TODO: What if bdhvn does not affect the amino acid. Check how this is handled
	private static List<String> determineSeqBDHVN(MutationSet<HIV> mutations) {
		MutationSet<HIV> bdhvnMuts = mutations.getAmbiguousCodons();
		return mutationListToTabularResult(bdhvnMuts);
	}

	private static List<String> determineSeqStopCodons(MutationSet<HIV> mutations) {
		MutationSet<HIV> stopCodons = mutations.getStopCodons();
		return mutationListToTabularResult(stopCodons);
	}

	private static List<String> determineSeqDeletions(MutationSet<HIV> mutations) {
		MutationSet<HIV> deletions = mutations.getDeletions();
		return mutationListToTabularResult(deletions);
	}

	private static List<String> determineSeqInsertions(MutationSet<HIV> mutations) {
		MutationSet<HIV> insertions = mutations.getInsertions();
		return mutationListToTabularResult(insertions);
	}


	private static List<String> determineFrameShiftText(AlignedSequence<HIV> alignedSeq) {
		List<FrameShift<HIV>> frameShifts = alignedSeq.getFrameShifts();
		List<String> frameShiftFields = new ArrayList<>();
		String frameShiftsString = FrameShift.joinFrameShifts(frameShifts);
		frameShiftFields.add(Integer.toString(frameShifts.size()));
		frameShiftFields.add(frameShiftsString);
		return frameShiftFields;
	}

	private static List<String> determineSdrms(AlignedSequence<HIV> alignedSeq) {
		List<String> sdrmList = new ArrayList<>();
		Map<Gene<HIV>, AlignedGeneSeq<HIV>> seqResult = alignedSeq.getAlignedGeneSequenceMap();
		Strain<HIV> strain = alignedSeq.getStrain();
		for (Gene<HIV> gene : strain.getGenes()) {
			String sdrmText = "NA";
			if (seqResult.containsKey(gene)) {
				MutationSet<HIV> sdrms = seqResult.get(gene).getSdrms();
				sdrmText = sdrms.join();
			}
			sdrmList.add(sdrmText);
		}
		return sdrmList;

	}

	public static List<String> determineMutLists(AlignedSequence<HIV> alignedSeq) {
		HIV hiv = HIV.getInstance();
		List<String> mutListStrings = new ArrayList<>();
		Map<Gene<HIV>, AlignedGeneSeq<HIV>> seqResult = alignedSeq.getAlignedGeneSequenceMap();
		Strain<HIV> strain = alignedSeq.getStrain();
		for (Gene<HIV> gene : strain.getGenes()) {
			if (!seqResult.containsKey(gene)) {
				if (gene.getAbstractGene().equals("RT")) {
					mutListStrings.addAll(Arrays.asList("NA", "NA", "NA", "NA", "NA"));
				}
				else {
					mutListStrings.addAll(Arrays.asList("NA", "NA", "NA"));
				}
			} else {
				AlignedGeneSeq<HIV> seq = seqResult.get(gene);
				List<MutationType<HIV>> mutTypes;
				if (gene.getName().equals("HIV1RT")) {
					mutTypes = Arrays.asList(
						hiv.getMutationType("NRTI"),
						null,
						hiv.getMutationType("NNRTI"),
						null,
						hiv.getMutationType("Other")
					);
				}
				else if (gene.getName().equals("HIV2ART") || gene.getName().equals("HIV2BRT")) {
					mutTypes = Arrays.asList(
						hiv.getMutationType("Major"),
						hiv.getMutationType("Accessory"),
						null,
						null,
						hiv.getMutationType("Other")
					);
				}
				else {
					mutTypes = Arrays.asList(
						hiv.getMutationType("Major"),
						hiv.getMutationType("Accessory"),
						hiv.getMutationType("Other")
					);
				}
				for (MutationType<HIV> mutType : mutTypes) {
					if (mutType == null) {
						mutListStrings.add("NA");
						continue;
					}
					MutationSet<HIV> mutTypeMutations = seq.getMutationsByMutType(mutType);
					if (mutTypeMutations.isEmpty()) {
						mutListStrings.add("None");
					} else {
						mutListStrings.add(mutTypeMutations.join());
					}
				}
			}
		}
		return mutListStrings;
	}


	// Four lists are returned. They should be in the following order: PI, NRTI, NNRTI, INSTI
	private static List<String> determineNonDrmTsms(AlignedSequence<HIV> alignedSeq) {
		HIV hiv = HIV.getInstance();
		List<String> nonDrmTsmsList = new ArrayList<>();
		Map<Gene<HIV>, AlignedGeneSeq<HIV>> seqResult = alignedSeq.getAlignedGeneSequenceMap();
		Strain<HIV> strain = alignedSeq.getStrain();

		for (Gene<HIV> gene : strain.getGenes()) {
			if (!seqResult.containsKey(gene)) {
				if (gene.getAbstractGene().equals("RT")) {
					nonDrmTsmsList.add("NA");
					nonDrmTsmsList.add("NA");
				} else {
					nonDrmTsmsList.add("NA");
				}
			} else {
				AlignedGeneSeq<HIV> seq = seqResult.get(gene);
				Map<DrugClass<HIV>, MutationSet<HIV>> allNonDrmTsms = seq.getNonDrmTsms();
				List<DrugClass<HIV>> drugClasses;
				switch (gene.getAbstractGene()) {
					case "PR":
						drugClasses = Arrays.asList(hiv.getDrugClass("PI"));
						break;
					case "RT":
						drugClasses = Arrays.asList(
							hiv.getDrugClass("NRTI"),
							hiv.getDrugClass("NNRTI")
						);
						break;
					default:
						// case IN:
						drugClasses = Arrays.asList(hiv.getDrugClass("INSTI"));
						break;
				}
				
				for (DrugClass<HIV> drugClass : drugClasses) {
					if (!allNonDrmTsms.containsKey(drugClass)) {
						nonDrmTsmsList.add("NA");
						continue;
					}
					MutationSet<HIV> nonDrmTsms = allNonDrmTsms.get(drugClass);

					String nonDrmTsmsText;
					if (nonDrmTsms.size() > 0) {
						nonDrmTsmsText = nonDrmTsms.join();
					} else {
						nonDrmTsmsText = "None";
					}
					nonDrmTsmsList.add(nonDrmTsmsText);
				}
			}
		}

		return nonDrmTsmsList;
	}


	private static List<String> determineGeneBoundaries (AlignedSequence<HIV> alignedSeq) {
		List<String>geneBoundaries = new ArrayList<>();
		Map<Gene<HIV>, AlignedGeneSeq<HIV>> seqResult = alignedSeq.getAlignedGeneSequenceMap();
		Strain<HIV> strain = alignedSeq.getStrain();
		for (Gene<HIV> gene : strain.getGenes()) {
			String firstAA = "NA", lastAA = "NA";
			if (seqResult.containsKey(gene)) {
				firstAA = "" + seqResult.get(gene).getFirstAA();
				lastAA = "" + seqResult.get(gene).getLastAA();
			}
			geneBoundaries.add(firstAA);
			geneBoundaries.add(lastAA);
		}
		return geneBoundaries;
	}

}
