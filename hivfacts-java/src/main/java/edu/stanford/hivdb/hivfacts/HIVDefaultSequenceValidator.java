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
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;


import com.google.common.collect.Lists;

import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.FrameShift;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.sequences.AlignedGeneSeq;
import edu.stanford.hivdb.sequences.AlignedSequence;
import edu.stanford.hivdb.sequences.SequenceValidator;
import edu.stanford.hivdb.utilities.Json;
import edu.stanford.hivdb.utilities.ValidationLevel;
import edu.stanford.hivdb.utilities.ValidationResult;

public class HIVDefaultSequenceValidator implements SequenceValidator<HIV> {

	private static final Map<String, ValidationLevel> VALIDATION_RESULT_LEVELS;
	private static final Map<String, String> VALIDATION_RESULT_MESSAGES;

	static {
		Map<String, ValidationLevel> levels = new HashMap<>();
		Map<String, String> messages = new HashMap<>();

		levels.put("no-gene-found", ValidationLevel.CRITICAL);
		messages.put("no-gene-found",
					"There were no Protease, Reverse Transcriptase, or " +
					"Integrase genes found, refuse to process.");

		levels.put("not-aligned-gene", ValidationLevel.SEVERE_WARNING);
		messages.put(
			"not-aligned-gene",
			"This sequence may also have nucleotides belonging to %s. " +
			"Analysis of this part of the input sequence was suppressed " +
			"due to poor quality, insufficient size or " +
			"improper concatenation of multiple partial sequences.");

		levels.put("gap-too-long", ValidationLevel.CRITICAL);
		messages.put(
			"gap-too-long",
			"This sequence has critical potentially correctable errors. It has a large sequence gap, " +
			"defined as an insertion or deletion of >30 bps. One possible cause of this error is that " +
			"the input sequence was concatenated from multiple partial sequences. Adding 'N's in place " +
			"of the missing sequence will allow the sequence to be processed.");

		levels.put("sequence-trimmed", ValidationLevel.WARNING);
		messages.put(
			"sequence-trimmed",
			"The %s sequence had %d amino acids trimmed from its %s-end due to poor quality.");

		levels.put("sequence-much-too-short", ValidationLevel.SEVERE_WARNING);
		messages.put(
			"sequence-much-too-short",
			"The %s sequence contains just %d codons, " +
			"which is not sufficient for a comprehensive interpretation.");

		levels.put("sequence-too-short", ValidationLevel.WARNING);
		messages.put(
			"sequence-too-short",
			"The %s sequence contains just %d codons, " +
			"which is not sufficient for a comprehensive interpretation.");

		levels.put("invalid-nas-removed", ValidationLevel.NOTE);
		messages.put(
			"invalid-nas-removed",
			"Non-NA character(s) %s were found and removed from the sequence.");

		levels.put("severe-warning-too-many-stop-codons", ValidationLevel.SEVERE_WARNING);
		messages.put("severe-warning-too-many-stop-codons", "There are %d stop codons in %s: %s.");

		levels.put("note-stop-codon", ValidationLevel.NOTE);
		messages.put("note-stop-codon", "There is %d stop codon in %s: %s.");

		levels.put("much-too-many-unusual-mutations", ValidationLevel.SEVERE_WARNING);
		messages.put("much-too-many-unusual-mutations", "There are %d unusual mutations in %s: %s.");

		levels.put("too-many-unusual-mutations", ValidationLevel.WARNING);
		messages.put("too-many-unusual-mutations", "There are %d unusual mutations in %s: %s.");

		levels.put("some-unusual-mutations", ValidationLevel.NOTE);
		messages.put("some-unusual-mutations", "There are %d unusual mutations in %s: %s.");

		levels.put("unusual-mutation-at-DRP-plural", ValidationLevel.WARNING);
		messages.put("unusual-mutation-at-DRP-plural",
				     "There are %d unusual mutations at drug-resistance positions in %s: %s.");

		levels.put("unusual-mutation-at-DRP", ValidationLevel.NOTE);
		messages.put("unusual-mutation-at-DRP",
				     "There is %d unusual mutation at a drug-resistance position in %s: %s.");

		levels.put("severe-APOBEC", ValidationLevel.SEVERE_WARNING);
		messages.put("severe-APOBEC", "The following %d APOBEC muts were present in the sequence: %s.%s");

		levels.put("definite-APOBEC", ValidationLevel.WARNING);
		messages.put("definite-APOBEC", "The following %d APOBEC muts were present in the sequence: %s.%s");

		levels.put("possible-APOBEC-influence", ValidationLevel.NOTE);
		messages.put("possible-APOBEC-influence", "The following %d APOBEC muts were present in the sequence: %s.%s");

		levels.put("multiple-apobec-at-DRP", ValidationLevel.SEVERE_WARNING);
		messages.put("multiple-apobec-at-DRP",
				     "There are %d APOBEC-associated mutations at drug-resistance positions: %s.");

		levels.put("single-apobec-at-DRP", ValidationLevel.WARNING);
		messages.put("single-apobec-at-DRP",
				     "There is %d APOBEC-associated mutation at a drug-resistance position: %s.");

		levels.put("two-or-more-unusual-indels-and-frameshifts", ValidationLevel.SEVERE_WARNING);
		messages.put("two-or-more-unusual-indels-and-frameshifts",
				"The %s gene has %d unusual indels and/or frameshifts. " +
				"The indels include %s. The frameshifts include %s.");

		levels.put("two-or-more-frameshifts", ValidationLevel.SEVERE_WARNING);
		messages.put("two-or-more-frameshifts", "The %s gene has %d frameshifts: %s.");

		levels.put("two-or-more-unusual-indels", ValidationLevel.SEVERE_WARNING);
		messages.put("two-or-more-unusual-indels", "The %s gene has %d unusual indels: %s.");

		levels.put("one-frameshift", ValidationLevel.WARNING);
		messages.put("one-frameshift", "The %s gene has a frameshift: %s.");

		levels.put("one-unusual-indel", ValidationLevel.WARNING);
		messages.put("one-unusual-indel", "The %s gene has an unusual indel: %s.");

		levels.put("hiv-2", ValidationLevel.WARNING);
		messages.put("hiv-2", "The sequence is from an HIV-2 virus");

		levels.put("overlap", ValidationLevel.WARNING);
		messages.put(
			"overlap", "Alignment overlap detected at the begining of %s " +
			"sequence (\"%s\"). Try insert Ns between partial sequences.");

		levels.put("reverse-complement", ValidationLevel.WARNING);
		messages.put(
			"reverse-complement", "This report was derived from the reverse complement of input sequence.");

		levels.put("unsequenced-region", ValidationLevel.WARNING);
		messages.put("unsequenced-region", "There are %d %s positions located in unsequenced region(s): %s.");

		VALIDATION_RESULT_LEVELS = Collections.unmodifiableMap(levels);
		VALIDATION_RESULT_MESSAGES = Collections.unmodifiableMap(messages);
	}
	
	protected HIVDefaultSequenceValidator() {}

	@Override
	public List<ValidationResult> validate(AlignedSequence<HIV> alignedSequence) {
		List<ValidationResult> results = new ArrayList<>();
		results.addAll(validateNotEmpty(alignedSequence));
		if (results.size() > 0) {
			return results;
		}
		results.addAll(validateReverseComplement(alignedSequence));
		results.addAll(validateGene(alignedSequence));
		results.addAll(validateSequenceSize(alignedSequence));
		results.addAll(validateUnsequencedRegion(alignedSequence));
		results.addAll(validateShrinkage(alignedSequence));
		results.addAll(validateLongGap(alignedSequence));
		results.addAll(validateNAs(alignedSequence));
		results.addAll(validateGaps(alignedSequence));
		results.addAll(validateNotApobec(alignedSequence));
		results.addAll(validateNoStopCodons(alignedSequence));
		results.addAll(validateNoTooManyUnusualMutations(alignedSequence));
		return results;
	}

	protected static ValidationResult newValidationResult(String key, Object... args) {
		ValidationLevel level = VALIDATION_RESULT_LEVELS.get(key);
		String message = String.format(
			VALIDATION_RESULT_MESSAGES.get(key),
			args);
		return new ValidationResult(level, message);
	}

	protected static List<ValidationResult> validateUnsequencedRegion(AlignedSequence<?> alignedSequence) {
		List<ValidationResult> results = new ArrayList<>();
		for (AlignedGeneSeq<?> geneSeq : alignedSequence.getAlignedGeneSequences()) {
			MutationSet<?> unsequenced = geneSeq.getMutations().filterBy(Mutation::isUnsequenced);
			if (unsequenced.size() > 2) {
				results.add(newValidationResult(
					"unsequenced-region", unsequenced.size(),
					geneSeq.getAbstractGene(), unsequenced.join(", ")
				));
			}
		}
		return results;
	}

	protected static List<ValidationResult> validateNotEmpty(AlignedSequence<?> alignedSequence) {
		if (alignedSequence.isEmpty()) {
			return Lists.newArrayList(newValidationResult("no-gene-found"));
		}
		return Collections.emptyList();
	}

	protected static List<ValidationResult> validateReverseComplement(AlignedSequence<?> alignedSequence) {
		if (alignedSequence.isReverseComplement()) {
			return Lists.newArrayList(newValidationResult("reverse-complement"));
		}
		return Collections.emptyList();
	}

	protected static List<ValidationResult> validateGene(AlignedSequence<?> alignedSequence) {
		Set<Gene<?>> discardedGenes = new LinkedHashSet<>(alignedSequence.getDiscardedGenes().keySet());
		int leftIgnored = 0x7fffffff;
		int rightIgnored = 0;
		Strain<?> strain = alignedSequence.getStrain();
		List<?> availableGenes = alignedSequence.getAvailableGenes();
		for (AlignedGeneSeq<?> geneSeq : alignedSequence.getAlignedGeneSequences()) {
			leftIgnored = Math.min(leftIgnored, geneSeq.getFirstNA() - 1);
			rightIgnored = Math.max(rightIgnored, geneSeq.getLastNA());
		}
		rightIgnored = alignedSequence.getInputSequence().getLength() - rightIgnored;
		if (!availableGenes.contains(strain.getGene("PR")) && leftIgnored > 210) {
			discardedGenes.add(strain.getGene("PR"));
		}
		if (!availableGenes.contains(strain.getGene("RT")) && leftIgnored > 800) {
			discardedGenes.add(strain.getGene("RT"));
		} else if (!availableGenes.contains(strain.getGene("RT")) && rightIgnored > 800) {
			discardedGenes.add(strain.getGene("RT"));
		}
		if (!availableGenes.contains(strain.getGene("IN")) && rightIgnored > 600) {
			discardedGenes.add(strain.getGene("IN"));
		}
		if (!discardedGenes.isEmpty()) {
			String textDiscardedGenes = discardedGenes
				.stream().map(g -> g.getName())
				.collect(Collectors.joining(" or "));
			return Lists.newArrayList(newValidationResult("not-aligned-gene", textDiscardedGenes));
		}
		return Collections.emptyList();
	}

	protected static List<ValidationResult> validateSequenceSize(AlignedSequence<HIV> alignedSequence) {
		int size;
		AlignedGeneSeq<?> geneSeq;
		int[] muchTooShortSize = new int[] {60, 150, 100};
		int[] tooShortSize = new int[] {80, 200, 200};
		String[] geneNames = new String[] {"PR", "RT", "IN"};
		List<ValidationResult> result = new ArrayList<>();
		
		for (int i = 0; i < 3; i ++) {
			geneSeq = alignedSequence.getAlignedGeneSequence(geneNames[i]);
			if (geneSeq != null) {
				size = geneSeq.getSize();
				if (size < muchTooShortSize[i]) {
					result.add(newValidationResult("sequence-much-too-short", geneNames[i], size));
				} else if (size < tooShortSize[i]) {
					result.add(newValidationResult("sequence-too-short", geneNames[i], size));
				}
			}
		}
		return result;
	}

	protected static List<ValidationResult> validateShrinkage(AlignedSequence<HIV> alignedSequence) {
		List<ValidationResult> result = new ArrayList<>();
		for (AlignedGeneSeq<HIV> geneSeq : alignedSequence.getAlignedGeneSequences()) {
			Gene<HIV> gene = geneSeq.getGene();
			int[] trimmed = geneSeq.getShrinkage();
			int leftTrimmed = trimmed[0];
			int rightTrimmed = trimmed[1];
			if (leftTrimmed > 0) {
				result.add(newValidationResult("sequence-trimmed", gene.getAbstractGene(), leftTrimmed, "5′"));
			}
			if (rightTrimmed > 0) {
				result.add(newValidationResult("sequence-trimmed", gene.getAbstractGene(), rightTrimmed, "3′"));
			}
		}
		return result;
	}

	protected static List<ValidationResult> validateLongGap(AlignedSequence<HIV> alignedSequence) {
		int gapLenThreshold = 10;
		int continuousDels = 0;
		List<ValidationResult> result = new ArrayList<>();
		for (Mutation<HIV> mut : alignedSequence.getMutations()) {
			if (continuousDels > gapLenThreshold) {
				result.add(newValidationResult("gap-too-long"));
				break;
			}
			if (mut.getInsertedNAs().length() > gapLenThreshold * 3) {
				result.add(newValidationResult("gap-too-long"));
				break;
			}
			if (mut.isDeletion()) {
				continuousDels ++;
			} else {
				continuousDels = 0;
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
			result.add(newValidationResult(
				"invalid-nas-removed", Json.dumps(String.join("", invalids)))
			);
		}
		return result;
	}

	protected static List<ValidationResult> validateNoStopCodons(AlignedSequence<HIV> alignedSequence) {
		Map<Gene<HIV>, AlignedGeneSeq<HIV>> alignedGeneSeqs = alignedSequence.getAlignedGeneSequenceMap();
		List<ValidationResult> result = new ArrayList<>();

		for (Gene<HIV> gene : alignedGeneSeqs.keySet()) {
			MutationSet<HIV> stopCodons = alignedGeneSeqs.get(gene).getStopCodons();
			String stops = stopCodons.join(", ");
			int numStopCodons = stopCodons.size();
			if (numStopCodons > 1) {
				result.add(newValidationResult(
					"severe-warning-too-many-stop-codons",
					numStopCodons, gene.getAbstractGene(), stops));
			} else if (numStopCodons > 0) {
				result.add(newValidationResult(
					"note-stop-codon", numStopCodons, gene.getAbstractGene(), stops));
			}
		}
		return result;
	}

	protected static List<ValidationResult> validateNoTooManyUnusualMutations(AlignedSequence<HIV> alignedSequence) {
		List<ValidationResult> result = new ArrayList<>();
		Map<Gene<HIV>, AlignedGeneSeq<HIV>> alignedGeneSeqs = alignedSequence.getAlignedGeneSequenceMap();

		for (Gene<HIV> gene : alignedGeneSeqs.keySet()) {
			AlignedGeneSeq<HIV> alignedGeneSeq = alignedGeneSeqs.get(gene);
			MutationSet<HIV> unusualMutations = alignedGeneSeq.getUnusualMutations();
			String text = unusualMutations.join(", ");
			int numUnusual = unusualMutations.size();
			if (numUnusual > 8) {
				result.add(newValidationResult(
					"much-too-many-unusual-mutations",
					numUnusual, gene.getAbstractGene(), text));
			} else if (numUnusual > 4) {
				result.add(newValidationResult(
					"too-many-unusual-mutations",
					numUnusual, gene.getAbstractGene(), text));
			} else if (numUnusual > 2) {
				result.add(newValidationResult(
					"some-unusual-mutations",
					numUnusual, gene.getAbstractGene(), text));
			}
			MutationSet<HIV> unusualMutAtDRP = alignedGeneSeq.getUnusualMutationsAtDrugResistancePositions();
			int numUnusualAtDRP = unusualMutAtDRP.size();
			if (numUnusualAtDRP > 1) {
				result.add(newValidationResult(
					"unusual-mutation-at-DRP-plural",
					numUnusualAtDRP, gene.getAbstractGene(), unusualMutAtDRP.join(", ")));
			} else if (numUnusualAtDRP == 1) {
				result.add(newValidationResult(
					"unusual-mutation-at-DRP",
					numUnusualAtDRP, gene.getAbstractGene(), unusualMutAtDRP.join(", ")));
			}
		}
		return result;
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

	protected static List<ValidationResult> validateNotApobec(AlignedSequence<HIV> alignedSequence) {
		MutationSet<HIV> apobecs = alignedSequence.getMutations().getApobecMutations();
		MutationSet<HIV> apobecDRMs = alignedSequence.getMutations().getApobecDRMs();
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
			results.add(newValidationResult("severe-APOBEC", numApobecMuts, apobecMutsText, extraCmt));
		} else if (numApobecMuts > 2) {
			results.add(newValidationResult("definite-APOBEC", numApobecMuts, apobecMutsText, extraCmt));
		} else if (numApobecMuts == 2) {
			results.add(newValidationResult("possible-APOBEC-influence", numApobecMuts, apobecMutsText, extraCmt));
		}

		MutationSet<HIV> apobecMutsAtDRP = apobecs.getAtDRPMutations();
		int numApobecMutsAtDRP = apobecMutsAtDRP.size();
		if (numApobecMutsAtDRP > 1) {
			results.add(newValidationResult(
				"multiple-apobec-at-DRP", numApobecMutsAtDRP,
				apobecMutsAtDRP.join(", ", Mutation::getHumanFormatWithGene)));
		} else if (numApobecMutsAtDRP == 1) {
			results.add(newValidationResult(
				"single-apobec-at-DRP", numApobecMutsAtDRP,
				apobecMutsAtDRP.join(", ", Mutation::getHumanFormatWithGene)));
		}
		return results;
	}

	private static List<ValidationResult> validateGaps(AlignedSequence<HIV> alignedSequence) {
		Map<Gene<HIV>, AlignedGeneSeq<HIV>> alignedGeneSeqs = alignedSequence.getAlignedGeneSequenceMap();
		List<Gene<HIV>> seqGenes = alignedSequence.getAvailableGenes();
		List<ValidationResult> results = new ArrayList<>();

		for (Gene<HIV> gene : seqGenes) {
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
					results.add(newValidationResult(
						"two-or-more-unusual-indels-and-frameshifts", gene.getAbstractGene(),
						numTotal, unusualIndelsListText, frameShiftListText));
				} else if (frameShifts.size() > 0) {
					results.add(newValidationResult(
						"two-or-more-frameshifts", gene.getAbstractGene(), numTotal,
						frameShiftListText));
				} else {
					results.add(newValidationResult(
						"two-or-more-unusual-indels", gene.getAbstractGene(), numTotal,
						unusualIndelsListText));
				}

			} else if (numTotal >0 ) {
				if (frameShifts.size() > 0) {
					results.add(newValidationResult(
						"one-frameshift", gene.getAbstractGene(), frameShiftListText));
				} else {
					results.add(newValidationResult(
						"one-unusual-indel", gene.getAbstractGene(), unusualIndelsListText));
				}

			}
		}
		return results;
	}



}
