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

package edu.stanford.hivdb.hivfacts.hiv2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import edu.stanford.hivdb.mutations.CodonReads;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.seqreads.GeneSequenceReads;
import edu.stanford.hivdb.seqreads.OneCodonReadsCoverage;
import edu.stanford.hivdb.seqreads.SequenceReads;
import edu.stanford.hivdb.seqreads.SequenceReadsValidator;
import edu.stanford.hivdb.utilities.ValidationLevel;
import edu.stanford.hivdb.utilities.ValidationResult;

public class HIV2DefaultSequenceReadsValidator implements SequenceReadsValidator<HIV2> {

	private static final Map<String, ValidationLevel> VALIDATION_RESULT_LEVELS;
	private static final Map<String, String> VALIDATION_RESULT_MESSAGES;
	// private static final Double PROPORTION_TRIMMED_POSITIONS_THRESHOLD = 0.05;
	// private static final Double GAP_LEN_THRESHOLD = 0.1;
	private static final Double UNUSUAL_THRESHOLD = 0.01;
	private static final Integer APOBEC_THRESHOLD = 2;
	// private static final Map<GeneEnum, Pair<Integer, Integer>> SEQUENCE_DRM_RANGES;
	// private static final Double SEQUENCE_DRM_MIN1 = 0.4;
	// private static final Double SEQUENCE_DRM_MIN2 = 0.6;

	static {
		Map<String, ValidationLevel> levels = new HashMap<>();
		Map<String, String> messages = new HashMap<>();

		levels.put("no-gene-found", ValidationLevel.CRITICAL);
		messages.put("no-gene-found",
					"There were no Protease, Reverse Transcriptase, or " +
					"Integrase genes found, refuse to process.");

		/* levels.put("gap-too-long", ValidationLevel.SEVERE_WARNING);
		messages.put(
			"gap-too-long",
			"More than 10%% of intermediate continuous positions are absent " +
			"in the submitted codon frequency table. This may indicates " +
			"preprocess issue and lead to incorrect subtyping detection. Click " +
			"the ‘Read Coverage’ button to review."); */
		
		levels.put("position-missing", ValidationLevel.WARNING);
		messages.put(
			"position-missing",
			"%d (%.1f%%) position%s %s missing between %s and %s. " +
			"Click the ‘Read Coverage’ button to review.");
		
		levels.put("drm-position-missing", ValidationLevel.WARNING);
		messages.put(
			"drm-position-missing",
			"%d (%.1f%%) out of %d position%s that may have a drug " +
			"resistance mutation %s missing between %s and %s. " +
			"Click the ‘Read Coverage’ button to review.");
		
		levels.put("min-read-depth-too-low", ValidationLevel.WARNING);
		messages.put(
			"min-read-depth-too-low",
			"You have selected a minimal read-depth of %d. However, " +
			"%d (%.1f%%) position%s in your sequence %s fewer than %d " +
			"reads. Click the ‘Read Coverage’ button to review.");

		/*levels.put("sequence-much-too-short", ValidationLevel.SEVERE_WARNING);
		messages.put(
			"sequence-much-too-short",
			"The %s sequence contains just %d codon(s) in DRM region, " +
			"which is not sufficient for a comprehensive interpretation.");

		levels.put("sequence-too-short", ValidationLevel.WARNING);
		messages.put(
			"sequence-too-short",
			"The %s sequence contains just %d codon(s) in DRM region, " +
			"which is not sufficient for a comprehensive interpretation.");*/

		levels.put("too-many-stop-codons", ValidationLevel.WARNING);
		messages.put("too-many-stop-codons", "There are %d stop codons: %s.");

		levels.put("one-stop-codon", ValidationLevel.WARNING);
		messages.put("one-stop-codon", "There is one stop codon: %s.");

		levels.put("too-many-unusual-mutations", ValidationLevel.SEVERE_WARNING);
		messages.put(
			"too-many-unusual-mutations",
			"At this threshold (%.1f%%), >= 1.0%% of positions have a " +
			"highly unusual mutation (defined as a prevalence <0.01%% in " +
			"published group M direct PCR sequences). This indicates that " +
			"there may be an unacceptably high risk that some mutations at " +
			"this threshold represent sequence artifacts.");

		
		levels.put("too-many-apobec-mutations-one-apobec-drm", ValidationLevel.WARNING);
		messages.put(
			"too-many-apobec-mutations-one-apobec-drm",
			"At this threshold (%.1f%%), >=3 positions with signature " +
			"APOBEC mutations. At this threshold, the sequence also contains " +
			"one drug-resistance mutation that could be caused by " +
			"APOBEC-mediated G-to-A hypermutation (%s). This " +
			"DRM therefore should be considered possible sequence artifacts.");
		
		levels.put("too-many-apobec-mutations-multiple-apobec-drms", ValidationLevel.WARNING);
		messages.put(
			"too-many-apobec-mutations-multiple-apobec-drms",
			"At this threshold (%.1f%%), >=3 positions with signature " +
			"APOBEC mutations. At this threshold, the sequence also contains " +
			"%d drug-resistance mutations that could be caused by " +
			"APOBEC-mediated G-to-A hypermutation (%s). These " +
			"DRMs therefore should be considered possible sequence artifacts.");

		levels.put("too-many-apobec-mutations-no-apobec-drm", ValidationLevel.WARNING);
		messages.put(
			"too-many-apobec-mutations-no-apobec-drm",
			"At this threshold (%.1f%%), >=3 positions with signature " +
			"APOBEC mutations. At this threshold, the sequence contains no" +
			"drug-resistance mutations that could be caused by " +
			"APOBEC-mediated G-to-A hypermutation.");

		levels.put("too-low-threshold", ValidationLevel.SEVERE_WARNING);
		messages.put(
			"too-low-threshold",
			"Extensive empirical data, as well as, modeling data suggest that " +
			"the risk of sequence artifact is almost inevitable at thresholds " +
			"below 0.5%% to 1.0%% unless unique molecular identifiers (UMI) are " +
			"used prior to PCR. We display quality control data for the 0.1%% " +
			"and 0.2%% thresholds solely to indicate that choosing a threshold " +
			"that is too low will result in sequence artifacts."
		);
		
		VALIDATION_RESULT_LEVELS = Collections.unmodifiableMap(levels);
		VALIDATION_RESULT_MESSAGES = Collections.unmodifiableMap(messages);

		/* Map<GeneEnum, Pair<Integer, Integer>> seqDrmRanges = new EnumMap<>(GeneEnum.class);
		seqDrmRanges.put(GeneEnum.PR, Pair.of(10, 90));
		seqDrmRanges.put(GeneEnum.RT, Pair.of(41, 348));
		seqDrmRanges.put(GeneEnum.IN, Pair.of(51, 263));
		SEQUENCE_DRM_RANGES = Collections.unmodifiableMap(seqDrmRanges); */

	}


	public List<ValidationResult> validate(SequenceReads<HIV2> seqReads) {
		List<ValidationResult> results = new ArrayList<>();
		results.addAll(validateNotEmpty(seqReads));
		if (!results.isEmpty()) {
			return results;
		}
		results.addAll(validateTooLowThreshold(seqReads));
		results.addAll(validateNoMissingPositions(seqReads));
		results.addAll(validateTrimmedPositions(seqReads));
		results.addAll(validateNoStopCodons(seqReads));
		results.addAll(validateNoTooManyUnusualMutations(seqReads));
		results.addAll(validateNoTooManyApobec(seqReads));
		return results;
	}

	protected ValidationResult newValidationResult(String key, Object... args) {
		ValidationLevel level = VALIDATION_RESULT_LEVELS.get(key);
		String message = String.format(
			VALIDATION_RESULT_MESSAGES.get(key),
			args);
		return new ValidationResult(level, message);
	}

	protected List<ValidationResult> validateNotEmpty(SequenceReads<HIV2> seqReads) {
		if (seqReads.isEmpty()) {
			return Lists.newArrayList(newValidationResult("no-gene-found"));
		}
		return Collections.emptyList();
	}
	
	protected List<ValidationResult> validateTooLowThreshold(SequenceReads<HIV2> seqReads) {
		List<ValidationResult> results = new ArrayList<>();
		double cutoff = seqReads.getMinPrevalence();
		if (Math.abs(cutoff - 0.001) < 1e-5 ||
			Math.abs(cutoff - 0.002) < 1e-5) {
			results.add(newValidationResult("too-low-threshold"));
		}
		return results;
				
	}
	
	protected List<ValidationResult> validateTrimmedPositions(SequenceReads<HIV2> seqReads) {
		List<ValidationResult> results = new ArrayList<>();
		List<OneCodonReadsCoverage<HIV2>> crcs = seqReads.getCodonReadsCoverage();
		long trimmedPos = crcs.stream().filter(crc -> crc.isTrimmed()).count();
		if (trimmedPos > 0) {
			double totalPos = crcs.size();
			double pcnt = (double) trimmedPos / totalPos;
			long minReadDepth = seqReads.getMinReadDepth();
			results.add(newValidationResult(
				"min-read-depth-too-low",
				minReadDepth, trimmedPos, pcnt * 100,
				trimmedPos == 1 ? "" : "s",
				trimmedPos == 1 ? "has" : "have",
				minReadDepth
			));
		}
		return results;
	}

	/*protected boolean validateSequenceSize() {
		long size;
		Gene[] genes = seqReads.getStrain().getGenes();
		GeneSequenceReads geneSeqReads;
		boolean validated = true;
		
		for (Gene gene : genes) {
			geneSeqReads = seqReads.getGeneSequenceReads(gene);
			if (geneSeqReads != null) {
				Pair<Integer, Integer> drmRange = SEQUENCE_DRM_RANGES.get(gene.getGeneEnum());
				Integer left = drmRange.getLeft();
				Integer right = drmRange.getRight();

				size = geneSeqReads.getAllPositionCodonReads().stream()
					.filter(pcr -> {
						long pos = pcr.getPosition();
						return pos >= left && pos <= right;
					})
					.count();
				double ratio = (double) size / (double) (right - left + 1);
				if (ratio < SEQUENCE_DRM_MIN1) {
					addValidationResult("sequence-much-too-short", gene, size);
					validated = false;
				}
				else if (ratio < SEQUENCE_DRM_MIN2) {
					addValidationResult("sequence-too-short", gene, size);
					validated = false;
				}
			}
		}
		return validated;
	}*/

	/*protected boolean validateLongGap() {
		boolean validated = true;
		List<OneCodonReadsCoverage> crcs = seqReads.getCodonReadsCoverage();
		if (crcs.isEmpty()) {
			return validated;
		}
		long leftMost = crcs.get(0).getPolPosition();
		long rightMost = crcs.get(crcs.size() - 1).getPolPosition();
		double maxGapWidth = Math.ceil((rightMost - leftMost + 1) * GAP_LEN_THRESHOLD);
		
		long prevPos = leftMost;
		for (OneCodonReadsCoverage crc : crcs) {
			long curPos = crc.getPolPosition();
			double posDiff = curPos - prevPos;
			if (posDiff > maxGapWidth) {
				addValidationResult("gap-too-long");
				validated = false;
				break;
			}
			prevPos = curPos;
		}
		for (Mutation mut : seqReads.getMutations()) {
			if (mut.getInsertedNAs().length() > maxGapWidth * 3) {
				addValidationResult("gap-too-long");
				validated = false;
				break;
			}
		}
		return validated;
	}*/
	
	protected List<ValidationResult> validateNoMissingPositions(SequenceReads<HIV2> seqReads) {
		List<ValidationResult> results = new ArrayList<>();
		List<OneCodonReadsCoverage<HIV2>> crcs = seqReads.getCodonReadsCoverage();
		if (crcs.isEmpty()) {
			return results;
		}
		GenePosition<HIV2> leftMost = crcs.get(0).getGenePosition();
		GenePosition<HIV2> rightMost = crcs.get(crcs.size() - 1).getGenePosition();
		Set<GenePosition<HIV2>> needGenePositions = (
			GenePosition.getGenePositionsBetween(leftMost, rightMost)
		);
		Set<GenePosition<HIV2>> needDRGenePositions = GenePosition.getDRGenePositionsBetween(leftMost, rightMost);
		Set<GenePosition<HIV2>> availableGenePositions = (
			crcs.stream().map(crc -> crc.getGenePosition()).collect(Collectors.toSet())
		);
		long totalPositions = needGenePositions.size();
		long totalDRPs = needDRGenePositions.size();
		long missingPositions = (
			needGenePositions.stream()
			.filter(gp -> !availableGenePositions.contains(gp))
			.count()
		);
		long missingDRPs = (
			needDRGenePositions.stream()
			.filter(gp -> !availableGenePositions.contains(gp))
			.count()
		);
		if (missingPositions > 0) {
			results.add(newValidationResult(
				"position-missing",
				missingPositions,
				(double) missingPositions / (double) totalPositions * 100.0,
				missingPositions == 1 ? "" : "s",
				missingPositions == 1 ? "is" : "are",
				leftMost.toString(), rightMost.toString()
			));
		}
		if (missingDRPs > 0) {
			results.add(newValidationResult(
				"drm-position-missing",
				missingDRPs,
				(double) missingDRPs / (double) totalDRPs * 100.0,
				totalDRPs,
				missingDRPs == 1 ? "" : "s",
				missingDRPs == 1 ? "is" : "are",
				leftMost.toString(), rightMost.toString()
			));
		}
		return results;
	}

	protected List<ValidationResult> validateNoStopCodons(SequenceReads<HIV2> seqReads) {
		List<ValidationResult> results = new ArrayList<>();
		List<GeneSequenceReads<HIV2>> allGeneSeqReads = seqReads.getAllGeneSequenceReads();
		MutationSet<HIV2> stopCodons = new MutationSet<>();

		for (GeneSequenceReads<HIV2> gsr : allGeneSeqReads) {
			stopCodons = stopCodons.mergesWith(gsr.getMutations().getStopCodons());
		}
		String stops = stopCodons.join(", ", Mutation::getHumanFormatWithGene);
		int numStopCodons = stopCodons.size();
		if (numStopCodons > 1) {
			results.add(newValidationResult(
				"too-many-stop-codons",
				numStopCodons, stops
			));
		} else if (numStopCodons > 0) {
			results.add(newValidationResult("one-stop-codon", stops));
		}
		return results;
	}

	protected List<ValidationResult> validateNoTooManyUnusualMutations(SequenceReads<HIV2> seqReads) {
		List<ValidationResult> results = new ArrayList<>();
		List<GeneSequenceReads<HIV2>> allGeneSeqReads = seqReads.getAllGeneSequenceReads();
		double numUnusuals = 0;
		double numPositions = 0;
		double cutoff = seqReads.getMinPrevalence();
		
		for (GeneSequenceReads<HIV2> gsr : allGeneSeqReads) {
			numUnusuals += (
				gsr.getAllPositionCodonReads()
				.stream()
				.mapToInt(pcr -> {
					for (CodonReads<HIV2> cr : pcr.getCodonReads(true, 1., cutoff)) {
						if (cr.isUnusual()) {
							return 1;
						}
					}
					return 0;
				})
				.sum());
			numPositions += gsr.getAllPositionCodonReads().size();
		}
		double unusualPcnt = numUnusuals / numPositions;
		if (unusualPcnt > UNUSUAL_THRESHOLD) {
			results.add(newValidationResult(
				"too-many-unusual-mutations", cutoff * 100
			));
		}
		return results;
	}

	protected List<ValidationResult> validateNoTooManyApobec(SequenceReads<HIV2> seqReads) {
		List<ValidationResult> results = new ArrayList<>();
		List<GeneSequenceReads<HIV2>> allGeneSeqReads = seqReads.getAllGeneSequenceReads();
		int numAPOBECs = 0;
		MutationSet<HIV2> apobecDRMs = new MutationSet<>();
		double cutoff = seqReads.getMinPrevalence();

		for (GeneSequenceReads<HIV2> gsr : allGeneSeqReads) {
			numAPOBECs += (
				gsr.getAllPositionCodonReads()
				.stream()
				.mapToInt(pcr -> {
					for (CodonReads<HIV2> cr : pcr.getCodonReads(true, 1., cutoff)) {
						if (cr.isApobecMutation()) {
							return 1;
						}
					}
					return 0;
				})
				.sum());
			apobecDRMs = apobecDRMs.mergesWith(gsr.getMutations().getApobecDRMs());
		}
		int numApobecDRMs = apobecDRMs.size();
		if (numAPOBECs > APOBEC_THRESHOLD) {
			String apobecs = apobecDRMs.join(", ", Mutation::getHumanFormatWithGene);
			
			if (numApobecDRMs > 1) {
				results.add(newValidationResult(
					"too-many-apobec-mutations-multiple-apobec-drms",
					cutoff * 100, numApobecDRMs, apobecs
				));
			}
			else if (numApobecDRMs == 1) {
				results.add(newValidationResult(
					"too-many-apobec-mutations-one-apobec-drm",
					cutoff * 100, apobecs
				));
			}
			else {
				results.add(newValidationResult(
					"too-many-apobec-mutations-no-apobec-drm",
					cutoff * 100
				));
			}
		}
		return results;
	}

}
