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
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;

import com.google.common.collect.Lists;

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.FrameShift;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.sequences.AlignedGeneSeq;
import edu.stanford.hivdb.sequences.AlignedSequence;
import edu.stanford.hivdb.sequences.SequenceValidator;
import edu.stanford.hivdb.sequences.GeneRegions;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.utilities.MyStringUtils;
import edu.stanford.hivdb.utilities.ValidationLevel;
import edu.stanford.hivdb.utilities.ValidationResult;

public class HIVDefaultSequenceValidator implements SequenceValidator<HIV> {

	protected HIVDefaultSequenceValidator() {}

	@Override
	public List<ValidationResult> validate(AlignedSequence<HIV> alignedSequence, Collection<String> includeGenes) {
		List<ValidationResult> results = new ArrayList<>();
		results.addAll(validateNotEmpty(alignedSequence, includeGenes));
		if (results.size() > 0) {
			return results;
		}
		results.addAll(validateReverseComplement(alignedSequence));
		results.addAll(validateNoMissingPositions(alignedSequence, includeGenes));
		results.addAll(validateShrinkage(alignedSequence, includeGenes));
		results.addAll(validateLongGap(alignedSequence, includeGenes));
		results.addAll(validateNAs(alignedSequence));
		results.addAll(validateGaps(alignedSequence, includeGenes));
		results.addAll(validateNotApobec(alignedSequence, includeGenes));
		results.addAll(validateNoStopCodons(alignedSequence, includeGenes));
		results.addAll(validateNoTooManyUnusualMutations(alignedSequence, includeGenes));
		return results;
	}

	/* protected static List<ValidationResult> validateUnsequencedRegion(
		AlignedSequence<?> alignedSequence,
		Collection<String> includeGenes
	) {
		List<ValidationResult> results = new ArrayList<>();
		for (AlignedGeneSeq<?> geneSeq : alignedSequence.getAlignedGeneSequences()) {
			String geneText = geneSeq.getAbstractGene();
			if (!includeGenes.contains(geneText)) {
				continue;
			}
			GeneRegions<?> unseqRegions = geneSeq.getUnsequencedRegions();
			long unseqSize = unseqRegions.getSize();
			if (unseqSize > 2) {
				results.add(HIV1ValidationMessage.FASTAUnsequencedRegion.format(
					unseqSize,
					geneText,
					MyStringUtils.andListFormat(unseqRegions.getRegions())
				));
			}
		}
		return results;
	} */

	protected static List<ValidationResult> validateNotEmpty(
		AlignedSequence<?> alignedSequence,
		Collection<String> includeGenes
	) {
		boolean isNotEmpty = alignedSequence.getAvailableGenes().stream()
			.anyMatch(gene -> includeGenes.contains(gene.getAbstractGene()));
		if (!isNotEmpty) {
			return Lists.newArrayList(HIV1ValidationMessage.NoGeneFound.format(
				MyStringUtils.andListFormat(includeGenes))
			);
		}
		return Collections.emptyList();
	}

	protected static List<ValidationResult> validateReverseComplement(AlignedSequence<?> alignedSequence) {
		if (alignedSequence.isReverseComplement()) {
			return Lists.newArrayList(HIV1ValidationMessage.FASTAReverseComplement.format());
		}
		return Collections.emptyList();
	}

	/* protected static List<ValidationResult> validateGene(
		AlignedSequence<?> alignedSequence,
		Collection<String> includeGenes
	) {
		Set<Gene<?>> discardedGenes = alignedSequence
			.getDiscardedGenes()
			.keySet()
			.stream()
			.filter(gene -> includeGenes.contains(gene.getAbstractGene()))
			.collect(Collectors.toCollection(LinkedHashSet::new));
		
		int leftIgnored = 0x7fffffff;
		int rightIgnored = 0;
		Strain<?> strain = alignedSequence.getStrain();
		List<?> availableGenes = alignedSequence.getAvailableGenes();
		for (AlignedGeneSeq<?> geneSeq : alignedSequence.getAlignedGeneSequences()) {
			leftIgnored = Math.min(leftIgnored, geneSeq.getFirstNA() - 1);
			rightIgnored = Math.max(rightIgnored, geneSeq.getLastNA());
		}
		rightIgnored = alignedSequence.getInputSequence().getLength() - rightIgnored;
		if (includeGenes.contains("PR")) {
			if (!availableGenes.contains(strain.getGene("PR")) && leftIgnored > 210) {
				discardedGenes.add(strain.getGene("PR"));
			}
		}
		if (includeGenes.contains("RT")) {
			if (!availableGenes.contains(strain.getGene("RT")) && leftIgnored > 800) {
				discardedGenes.add(strain.getGene("RT"));
			} else if (!availableGenes.contains(strain.getGene("RT")) && rightIgnored > 800) {
				discardedGenes.add(strain.getGene("RT"));
			}
		}
		if (includeGenes.contains("IN")) {
			if (!availableGenes.contains(strain.getGene("IN")) && rightIgnored > 600) {
				discardedGenes.add(strain.getGene("IN"));
			}
		}
		if (!discardedGenes.isEmpty()) {
			String textDiscardedGenes = MyStringUtils.orListFormat(
				discardedGenes
					.stream().map(g -> g.getAbstractGene())
					.collect(Collectors.toList())
			);
			return Lists.newArrayList(HIV1ValidationMessage.FASTAGeneNotAligned.format(textDiscardedGenes));
		}
		return Collections.emptyList();
	} */

	/* protected static List<ValidationResult> validateSequenceSize(
		AlignedSequence<HIV> alignedSequence,
		Collection<String> includeGenes
	) {
		int size;
		AlignedGeneSeq<?> geneSeq;
		int[] muchTooShortSize = new int[] {60, 150, 100};
		int[] tooShortSize = new int[] {80, 200, 200};
		String[] geneNames = new String[] {"PR", "RT", "IN"};
		List<ValidationResult> result = new ArrayList<>();
		
		for (int i = 0; i < 3; i ++) {
			String geneName = geneNames[i];
			if (!includeGenes.contains(geneName)) {
				continue;
			}
			geneSeq = alignedSequence.getAlignedGeneSequence(geneNames[i]);
			if (geneSeq != null) {
				size = geneSeq.getSize();
				if (size < muchTooShortSize[i]) {
					result.add(HIV1ValidationMessage.FASTASequenceTooShort.formatWithLevel(
						ValidationLevel.SEVERE_WARNING,
						geneNames[i],
						size
					));
				} else if (size < tooShortSize[i]) {
					result.add(HIV1ValidationMessage.FASTASequenceTooShort.formatWithLevel(
						ValidationLevel.WARNING,	
						geneNames[i],
						size
					));
				}
			}
		}
		return result;
	} */
	
	protected static List<ValidationResult> validateNoMissingPositions(
		final Set<GenePosition<HIV>> needGenePositions,
		final Set<GenePosition<HIV>> needDRGenePositions,
		final Set<GenePosition<HIV>> availableGenePositions
	) {
		List<ValidationResult> results = new ArrayList<>();

		List<GenePosition<HIV>> missingPositions = needGenePositions.stream()
				.filter(gp -> !availableGenePositions.contains(gp))
				.collect(Collectors.toList());
		long numMissingPositions = missingPositions.size();

		List<GenePosition<HIV>> missingDRPs = needDRGenePositions.stream()
				.filter(gp -> !availableGenePositions.contains(gp))
				.collect(Collectors.toList());
		long numMissingDRPs = missingDRPs.size();
		
		String textMissingPositions = StringUtils.join(
			GeneRegions.newListOfGeneRegions(missingPositions),
			"; "
		);

		String textMissingDRPs = StringUtils.join(
			GeneRegions.newListOfGeneRegions(missingDRPs),
			"; "
		);
		
		if (numMissingDRPs > 1) {
			results.add(HIV1ValidationMessage.MultiplePositionsMissingWithMultipleDRPs.formatWithLevel(
				numMissingDRPs > 5 ? ValidationLevel.SEVERE_WARNING : (
					numMissingDRPs > 3 ? ValidationLevel.WARNING : ValidationLevel.NOTE
				),
				numMissingPositions,
				textMissingPositions,
				numMissingDRPs,
				textMissingDRPs
			));
		}
		else if (numMissingDRPs > 0 && numMissingPositions > 1) {
			results.add(HIV1ValidationMessage.MultiplePositionsMissingWithSingleDRP.formatWithLevel(
				ValidationLevel.WARNING,
				numMissingPositions,
				textMissingPositions,
				textMissingDRPs
			));
		}
		else if (numMissingDRPs > 0) {
			results.add(HIV1ValidationMessage.SingleDRPMissing.formatWithLevel(
				ValidationLevel.NOTE,
				textMissingDRPs
			));
		}
		else if (numMissingPositions > 1) {
			results.add(HIV1ValidationMessage.MultiplePositionsMissingWithoutDRP.formatWithLevel(
				ValidationLevel.WARNING,
				numMissingPositions,
				textMissingPositions
			));
		}
		else if (numMissingPositions > 0) {
			results.add(HIV1ValidationMessage.SinglePositionMissingWithoutDRP.formatWithLevel(
				ValidationLevel.NOTE,
				textMissingPositions
			));
		}
		return results;
	}

	protected static List<ValidationResult> validateNoMissingPositions(
			AlignedSequence<HIV> alignedSequence,
		Collection<String> includeGenes
	) {
		List<AlignedGeneSeq<HIV>> geneSeqs = alignedSequence.getAlignedGeneSequences(includeGenes);
		if (geneSeqs.isEmpty()) {
			return Collections.emptyList();
		}
		AlignedGeneSeq<HIV> geneSeq = geneSeqs.get(0);
		GenePosition<HIV> leftMost = new GenePosition<>(geneSeq.getGene(), 1);
		geneSeq = geneSeqs.get(geneSeqs.size() - 1);
		GenePosition<HIV> rightMost = new GenePosition<>(geneSeq.getGene(), geneSeq.getGene().getAASize());
		Strain<HIV> strain = alignedSequence.getStrain();
		Map<Gene<HIV>, GeneRegions<HIV>> unseqRegions = includeGenes.stream()
			.map(absGene -> strain.getGene(absGene))
			.collect(Collectors.toMap(
				gene -> gene,
				gene -> {
					AlignedGeneSeq<HIV> gs = alignedSequence.getAlignedGeneSequence(gene);
					return gs == null ? (
						GeneRegions.newGeneRegions(gene, 1, gene.getAASize())
					) : gs.getUnsequencedRegions();
				}
			));
			
		Set<GenePosition<HIV>> needGenePositions = GenePosition
			.getGenePositionsBetween(leftMost, rightMost, includeGenes);
		
		// For DRPs, the leftMost must be the begining of the first gene and the rightMost must be the ending of the last gene
		Set<GenePosition<HIV>> needDRGenePositions = GenePosition
			.getDRGenePositionsBetween(leftMost, rightMost, includeGenes);

		Set<GenePosition<HIV>> availableGenePositions = needGenePositions.stream()
			.filter(gpos -> {
				Gene<HIV> gene = gpos.getGene();
				GeneRegions<HIV> geneUnseqRegions = unseqRegions.get(gene);
				if (geneUnseqRegions == null) {
					return true;
				}
				if (geneUnseqRegions.contains(gpos.getPosition())) {
					return false;
				}
				return true;
			})	
			.collect(Collectors.toSet());
		return validateNoMissingPositions(
			needGenePositions,
			needDRGenePositions,
			availableGenePositions
		);
	}

	protected static List<ValidationResult> validateShrinkage(
		AlignedSequence<HIV> alignedSequence,
		Collection<String> includeGenes
	) {
		List<ValidationResult> result = new ArrayList<>();
		for (AlignedGeneSeq<HIV> geneSeq : alignedSequence.getAlignedGeneSequences()) {
			String geneText = geneSeq.getAbstractGene();
			if (!includeGenes.contains(geneText)) {
				continue;
			}
			int[] trimmed = geneSeq.getShrinkage();
			int leftTrimmed = trimmed[0];
			int rightTrimmed = trimmed[1];
			if (leftTrimmed > 0) {
				result.add(HIV1ValidationMessage.FASTASequenceTrimmed.format(geneText, leftTrimmed, "5′"));
			}
			if (rightTrimmed > 0) {
				result.add(HIV1ValidationMessage.FASTASequenceTrimmed.format(geneText, rightTrimmed, "3′"));
			}
		}
		return result;
	}

	protected static List<ValidationResult> validateLongGap(
		AlignedSequence<HIV> alignedSequence,
		Collection<String> includeGenes
	) {
		int gapLenThreshold = 10;
		int totalIndels = 0;
		List<ValidationResult> result = new ArrayList<>();
		for (Mutation<HIV> mut : alignedSequence.getMutations()) {
			if (!includeGenes.contains(mut.getAbstractGene())) {
				continue;
			}
			if (totalIndels > gapLenThreshold) {
				result.add(HIV1ValidationMessage.FASTAGapTooLong.format());
				break;
			}
			if (mut.getInsertedNAs().length() > gapLenThreshold * 3) {
				result.add(HIV1ValidationMessage.FASTAGapTooLong.format());
				break;
			}
			if (mut.isDeletion()) {
				totalIndels ++;
			}
			else if (mut.isInsertion()) {
				totalIndels += Math.round(mut.getInsertedNAs().length() / 3);
			}
			else {
				totalIndels = 0;
			}
		}
		return result;
	}

	protected static List<ValidationResult> validateNAs(AlignedSequence<HIV> alignedSequence) {
		List<String> invalids =
			alignedSequence.getInputSequence().removedInvalidChars()
			.stream().map(c -> "" + c)
			.collect(Collectors.toList());
		List<ValidationResult> result = new ArrayList<>();
		if (!invalids.isEmpty()) {
			result.add(HIV1ValidationMessage.FASTAInvalidNAsRemoved.format(
				Json.dumps(String.join("", invalids))
			));
		}
		return result;
	}

	protected static List<ValidationResult> validateNoStopCodons(
		AlignedSequence<HIV> alignedSequence,
		Collection<String> includeGenes
	) {
		List<ValidationResult> results = new ArrayList<>();
		MutationSet<HIV> stopCodons = (
			alignedSequence.getMutations()
			.getStopCodons()
			.filterBy(mut -> includeGenes.contains(mut.getAbstractGene()))
		);
		for (Map.Entry<Gene<HIV>, MutationSet<HIV>> entry : stopCodons.groupByGene().entrySet()) {
			String geneText = entry.getKey().getAbstractGene();
			MutationSet<HIV> geneStopCodons = entry.getValue();
			int numGeneStopCodons = geneStopCodons.size();
			String geneStopText = geneStopCodons.join(", ", Mutation::getHumanFormatWithAbstractGene);
			if (numGeneStopCodons > 1) {
				results.add(HIV1ValidationMessage.MultipleStopCodons.formatWithLevel(
					ValidationLevel.SEVERE_WARNING,
					numGeneStopCodons,
					geneText,
					geneStopText
				));
			} else if (numGeneStopCodons > 0) {
				results.add(HIV1ValidationMessage.SingleStopCodon.formatWithLevel(
					ValidationLevel.NOTE,
					geneText,
					geneStopText
				));
			}
		}
		
		return results;
	}

	protected static List<ValidationResult> validateNoTooManyUnusualMutations(
		AlignedSequence<HIV> alignedSequence,
		Collection<String> includeGenes
	) {
		List<ValidationResult> results = new ArrayList<>();
		MutationSet<HIV> unusualMuts = alignedSequence
			.getMutations()
			.getUnusualMutations()
			.filterBy(mut -> includeGenes.contains(mut.getAbstractGene()));
		
		for (Map.Entry<Gene<HIV>, MutationSet<HIV>> entry : unusualMuts.groupByGene().entrySet()) {
			String geneText = entry.getKey().getAbstractGene();
			MutationSet<HIV> geneUnusualMuts = entry.getValue();
			String geneMutText = geneUnusualMuts.join(", ");
			int numGeneUnusual = geneUnusualMuts.size();
			if (numGeneUnusual > 8) {
				results.add(HIV1ValidationMessage.MultipleUnusualMutations.formatWithLevel(
					ValidationLevel.SEVERE_WARNING,
					numGeneUnusual,
					geneText,
					geneMutText
				));
			} else if (numGeneUnusual > 4) {
				results.add(HIV1ValidationMessage.MultipleUnusualMutations.formatWithLevel(
					ValidationLevel.WARNING,
					numGeneUnusual,
					geneText,
					geneMutText
				));
			} else if (numGeneUnusual > 2) {
				results.add(HIV1ValidationMessage.MultipleUnusualMutations.formatWithLevel(
					ValidationLevel.NOTE,
					numGeneUnusual,
					geneText,
					geneMutText
				));
			}
			MutationSet<HIV> geneUnusualMutsAtDRP = geneUnusualMuts.getAtDRPMutations();
			int numGeneUnusualAtDRP = geneUnusualMutsAtDRP.size();
			if (numGeneUnusualAtDRP > 1) {
				results.add(HIV1ValidationMessage.MultipleUnusualMutationsAtDRP.formatWithLevel(
					ValidationLevel.WARNING,
					numGeneUnusualAtDRP,
					geneText,
					geneUnusualMutsAtDRP.join(", ")
				));
			} else if (numGeneUnusualAtDRP == 1) {
				results.add(HIV1ValidationMessage.SingleUnusualMutationAtDRP.formatWithLevel(
					ValidationLevel.NOTE,
					geneText,
					geneUnusualMutsAtDRP.join(", ")
				));
			}
		}
		return results;
	}

	/*@SuppressWarnings("unused")
	private boolean validateGroupO() {
		Subtyper subtyper = alignedSequence.getSubtyper();
		return subtyper.getClosestSubtype() == Subtype.O;
	}

	@SuppressWarnings("unused")
	private boolean validateGroupN() {
		Subtyper subtyper = alignedSequence.getSubtyper();
		return subtyper.getClosestSubtype() == Subtype.N;
	}*/

	protected static List<ValidationResult> validateNotApobec(
		AlignedSequence<HIV> alignedSequence,
		Collection<String> includeGenes
	) {
		MutationSet<HIV> apobecs = alignedSequence
			.getMutations()
			.getApobecMutations()
			.filterByNoSplit(mut -> includeGenes.contains(mut.getAbstractGene()));
		MutationSet<HIV> apobecDRMs = alignedSequence
			.getMutations()
			.getApobecDRMs()
			.filterByNoSplit(mut -> includeGenes.contains(mut.getAbstractGene()));
		List<ValidationResult> results = new ArrayList<>();
		int numApobecMuts = apobecs.size();
		int numApobecDRMs = apobecDRMs.size();
		String apobecMutsText = (
			apobecs.groupByGene()
			.entrySet()
			.stream()
			.map(e -> String.format(
				"%s: %s", e.getKey().getAbstractGene(),
				e.getValue().join(", ")
			))
			.collect(Collectors.joining("; "))
		);
		String extraCmt = "";
		
		if (numApobecDRMs > 0) {
			extraCmt = String.format(
				" The following %d DRMs in this sequence could reflect APOBEC activity: %s.",
				numApobecDRMs,
				apobecDRMs.groupByGene()
				.entrySet()
				.stream()
				.map(e -> String.format(
					"%s: %s", e.getKey().getAbstractGene(),
					e.getValue().join(", ")
				))
				.collect(Collectors.joining("; "))
			);
		}

		if (numApobecMuts > 4) {
			results.add(HIV1ValidationMessage.MultipleApobec.formatWithLevel(
				ValidationLevel.SEVERE_WARNING,
				numApobecMuts,
				apobecMutsText,
				extraCmt
			));
		} else if (numApobecMuts > 2) {
			results.add(HIV1ValidationMessage.MultipleApobec.formatWithLevel(
				ValidationLevel.WARNING,	
				numApobecMuts,
				apobecMutsText,
				extraCmt
			));
		} else if (numApobecMuts == 2) {
			results.add(HIV1ValidationMessage.SingleApobec.formatWithLevel(
				ValidationLevel.NOTE,
				apobecMutsText,
				extraCmt
			));
		}

		MutationSet<HIV> apobecMutsAtDRP = apobecs.getAtDRPMutations();
		int numApobecMutsAtDRP = apobecMutsAtDRP.size();
		if (numApobecMutsAtDRP > 1) {
			results.add(HIV1ValidationMessage.MultipleApobecAtDRP.formatWithLevel(
				ValidationLevel.SEVERE_WARNING,
				numApobecMutsAtDRP,
				apobecMutsAtDRP.join(", ", Mutation::getHumanFormatWithAbstractGene)
			));
		} else if (numApobecMutsAtDRP == 1) {
			results.add(HIV1ValidationMessage.SingleApobecAtDRP.formatWithLevel(
				ValidationLevel.WARNING,
				apobecMutsAtDRP.join(", ", Mutation::getHumanFormatWithAbstractGene)
			));
		}
		return results;
	}

	private static List<ValidationResult> validateGaps(
		AlignedSequence<HIV> alignedSequence,
		Collection<String> includeGenes
	) {
		Map<Gene<HIV>, AlignedGeneSeq<HIV>> alignedGeneSeqs = alignedSequence.getAlignedGeneSequenceMap();
		List<Gene<HIV>> seqGenes = alignedSequence.getAvailableGenes();
		List<ValidationResult> results = new ArrayList<>();

		for (Gene<HIV> gene : seqGenes) {
			String geneText = gene.getAbstractGene();
			if (!includeGenes.contains(geneText)) {
				continue;
			}
			AlignedGeneSeq<HIV> alignedGeneSeq = alignedGeneSeqs.get(gene);
			List<FrameShift<HIV>> frameShifts = alignedGeneSeq.getFrameShifts();
			MutationSet<HIV> insertions = alignedGeneSeq.getInsertions();
			MutationSet<HIV> deletions = alignedGeneSeq.getDeletions();
			MutationSet<HIV> unusualInsertions = insertions.getUnusualMutations();
			MutationSet<HIV> unusualDeletions = deletions.getUnusualMutations();
			MutationSet<HIV> unusualIndels = unusualInsertions.mergesWith(unusualDeletions);
			int numTotal = frameShifts.size() + unusualInsertions.size() + unusualDeletions.size();
			String frameShiftListText = FrameShift.joinFrameShifts(frameShifts);
			String unusualIndelsListText = unusualIndels.join(", ");

			if (numTotal > 1) {
				if (frameShifts.size() > 0 && unusualIndels.size() > 0) {
					results.add(HIV1ValidationMessage.MultipleUnusualIndelsAndFrameshifts.formatWithLevel(
						ValidationLevel.SEVERE_WARNING,
						geneText,
						numTotal,
						unusualIndelsListText,
						frameShiftListText
					));
				} else if (frameShifts.size() > 0) {
					results.add(HIV1ValidationMessage.MultipleFrameShifts.formatWithLevel(
						ValidationLevel.SEVERE_WARNING,
						geneText,
						numTotal,
						frameShiftListText
					));
				} else {
					results.add(HIV1ValidationMessage.MultipleUnusualIndels.formatWithLevel(
						ValidationLevel.SEVERE_WARNING,
						geneText,
						numTotal,
						unusualIndelsListText
					));
				}

			} else if (numTotal >0 ) {
				if (frameShifts.size() > 0) {
					results.add(HIV1ValidationMessage.SingleFrameshift.formatWithLevel(
						ValidationLevel.WARNING,
						geneText,
						frameShiftListText
					));
				} else {
					results.add(HIV1ValidationMessage.SingleUnusualIndel.formatWithLevel(
						ValidationLevel.WARNING,
						geneText,
						unusualIndelsListText
					));
				}

			}
		}
		return results;
	}

}
