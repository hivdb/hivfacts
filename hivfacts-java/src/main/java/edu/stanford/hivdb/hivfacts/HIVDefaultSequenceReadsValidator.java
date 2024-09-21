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

import com.google.common.collect.Lists;

import edu.stanford.hivdb.mutations.CodonReads;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.GenePosition;
import edu.stanford.hivdb.seqreads.GeneSequenceReads;
import edu.stanford.hivdb.seqreads.OneCodonReadsCoverage;
import edu.stanford.hivdb.seqreads.SequenceReads;
import edu.stanford.hivdb.seqreads.SequenceReadsValidator;
import edu.stanford.hivdb.sequences.GeneRegions;
import edu.stanford.hivdb.utilities.MyStringUtils;
import edu.stanford.hivdb.utilities.ValidationLevel;
import edu.stanford.hivdb.utilities.ValidationResult;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;

public class HIVDefaultSequenceReadsValidator implements SequenceReadsValidator<HIV> {

	// private static final Double PROPORTION_TRIMMED_POSITIONS_THRESHOLD = 0.05;
	// private static final Double GAP_LEN_THRESHOLD = 0.1;
	private static final Double UNUSUAL_THRESHOLD = 0.01;
	private static final Integer APOBEC_THRESHOLD = 2;
	// private static final Map<GeneEnum, Pair<Integer, Integer>> SEQUENCE_DRM_RANGES;
	// private static final Double SEQUENCE_DRM_MIN1 = 0.4;
	// private static final Double SEQUENCE_DRM_MIN2 = 0.6;
	
	public List<ValidationResult> validate(SequenceReads<HIV> seqReads, Collection<String> includeGenes) {
		List<ValidationResult> results = new ArrayList<>();
		results.addAll(validateNotEmpty(seqReads, includeGenes));
		if (!results.isEmpty()) {
			return results;
		}
		results.addAll(validateTooLowThreshold(seqReads));
		results.addAll(validateNoMissingPositions(seqReads, includeGenes));
		results.addAll(validateTrimmedPositions(seqReads, includeGenes));
		results.addAll(validateNoStopCodons(seqReads, includeGenes));
		results.addAll(validateNoTooManyUnusualMutations(seqReads, includeGenes));
		results.addAll(validateNoTooManyApobec(seqReads, includeGenes));
		return results;
	}

	protected static List<ValidationResult> validateNotEmpty(
		SequenceReads<HIV> seqReads,
		Collection<String> includeGenes
	) {
		boolean isNotEmpty = seqReads.getAvailableGenes().stream()
			.anyMatch(gene -> includeGenes.contains(gene.getAbstractGene()));
		if (!isNotEmpty) {
			return Lists.newArrayList(
				HIV1ValidationMessage.NoGeneFound.format(MyStringUtils.andListFormat(includeGenes))
			);
		}
		return Collections.emptyList();
	}
	
	protected static List<ValidationResult> validateTooLowThreshold(SequenceReads<HIV> seqReads) {
		List<ValidationResult> results = new ArrayList<>();
		double cutoff = seqReads.getMinPrevalence();
		if (Math.abs(cutoff - 0.001) < 1e-5 ||
			Math.abs(cutoff - 0.002) < 1e-5) {
			results.add(HIV1ValidationMessage.NGSTooLowThreshold.format());
		}
		return results;
				
	}
	
	protected static List<ValidationResult> validateTrimmedPositions(
		SequenceReads<HIV> seqReads,
		Collection<String> includeGenes
	) {
		List<ValidationResult> results = new ArrayList<>();
		List<OneCodonReadsCoverage<HIV>> crcs = seqReads.getCodonReadsCoverage(includeGenes);
		long trimmedPos = crcs.stream().filter(crc -> crc.isTrimmed()).count();
		if (trimmedPos > 0) {
			double totalPos = crcs.size();
			double pcnt = (double) trimmedPos / totalPos;
			long minReadDepth = seqReads.getMinPositionReads();
			results.add(HIV1ValidationMessage.NGSMinReadDepthTooLow.format(
				minReadDepth, trimmedPos, pcnt * 100,
				trimmedPos == 1 ? "" : "s",
				trimmedPos == 1 ? "has" : "have",
				minReadDepth
			));
		}
		return results;
	}

	protected static List<ValidationResult> validateNoMissingPositions(
		SequenceReads<HIV> seqReads,
		Collection<String> includeGenes
	) {
		List<ValidationResult> results = new ArrayList<>();
		List<OneCodonReadsCoverage<HIV>> crcs = seqReads.getCodonReadsCoverage(includeGenes);
		if (crcs.isEmpty()) {
			return results;
		}
		OneCodonReadsCoverage<HIV> tmpCRC = crcs.get(0);
		GenePosition<HIV> leftMost = new GenePosition<>(tmpCRC.getGene(), 1);
		tmpCRC = crcs.get(crcs.size() - 1);
		GenePosition<HIV> rightMost = new GenePosition<>(tmpCRC.getGene(), tmpCRC.getGene().getAASize());
		
		Set<GenePosition<HIV>> needGenePositions = GenePosition
			.getGenePositionsBetween(leftMost, rightMost, includeGenes);

		// For DRPs, the leftMost must be the begining of the first gene and the rightMost must be the ending of the last gene
		Set<GenePosition<HIV>> needDRGenePositions = GenePosition
			.getDRGenePositionsBetween(leftMost, rightMost, includeGenes);
	
		Strain<HIV> strain = seqReads.getStrain();
		Map<Gene<HIV>, GeneRegions<HIV>> unseqRegions = includeGenes.stream()
			.map(absGene -> strain.getGene(absGene))
			.collect(Collectors.toMap(
				gene -> gene,
				gene -> {
					GeneSequenceReads<HIV> gs = seqReads.getGeneSequenceReads(gene);
					return gs == null ? (
						GeneRegions.newGeneRegions(gene, 1, gene.getAASize())
					) : gs.getUnsequencedRegions();
				}
			));

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
		
		return HIVDefaultSequenceValidator.validateNoMissingPositions(
				needGenePositions,
				needDRGenePositions,
				availableGenePositions
		);
	}

	protected static List<ValidationResult> validateNoStopCodons(
		SequenceReads<HIV> seqReads,
		Collection<String> includeGenes
	) {
		List<ValidationResult> results = new ArrayList<>();
		MutationSet<HIV> stopCodons = (
			seqReads.getMutations()
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
					ValidationLevel.WARNING,
					numGeneStopCodons,
					geneText,
					geneStopText
				));
			} else if (numGeneStopCodons > 0) {
				results.add(HIV1ValidationMessage.SingleStopCodon.formatWithLevel(
					ValidationLevel.WARNING,
					geneText,
					geneStopText
				));
			}
		}
		
		return results;
	}

	protected static List<ValidationResult> validateNoTooManyUnusualMutations(
		SequenceReads<HIV> seqReads,
		Collection<String> includeGenes
	) {
		List<ValidationResult> results = new ArrayList<>();
		List<GeneSequenceReads<HIV>> allGeneSeqReads = seqReads.getAllGeneSequenceReads(includeGenes);
		double numUnusuals = 0;
		double numPositions = 0;
		double cutoff = seqReads.getMinPrevalence();
		
		for (GeneSequenceReads<HIV> gsr : allGeneSeqReads) {
			numUnusuals += (
				gsr.getAllPositionCodonReads()
				.stream()
				.mapToInt(pcr -> {
					for (CodonReads<HIV> cr : pcr.getCodonReads(true, 1., cutoff)) {
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
			results.add(HIV1ValidationMessage.NGSTooManyUnusualMutations.format(cutoff * 100));
		}
		return results;
	}

	protected static List<ValidationResult> validateNoTooManyApobec(
		SequenceReads<HIV> seqReads,
		Collection<String> includeGenes
	) {
		List<ValidationResult> results = new ArrayList<>();
		List<GeneSequenceReads<HIV>> allGeneSeqReads = seqReads.getAllGeneSequenceReads(includeGenes);
		int numAPOBECs = 0;
		MutationSet<HIV> apobecDRMs = new MutationSet<>();
		double cutoff = seqReads.getMinPrevalence();

		for (GeneSequenceReads<HIV> gsr : allGeneSeqReads) {
			numAPOBECs += (
				gsr.getAllPositionCodonReads()
				.stream()
				.mapToInt(pcr -> {
					for (CodonReads<HIV> cr : pcr.getCodonReads(true, 1., cutoff)) {
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
			String apobecs = apobecDRMs.join(", ", Mutation::getHumanFormatWithAbstractGene);
			
			if (numApobecDRMs > 1) {
				results.add(HIV1ValidationMessage.NGSTooManyApobecMutationsMultipleApobecDRMs.format(
					cutoff * 100, numApobecDRMs, apobecs
				));
			}
			else if (numApobecDRMs == 1) {
				results.add(HIV1ValidationMessage.NGSTooManyApobecMutationsOneApobecDRM.format(
					cutoff * 100, apobecs
				));
			}
			else {
				results.add(HIV1ValidationMessage.NGSTooManyApobecMutationsNoApobecDRM.format(
					cutoff * 100
				));
			}
		}
		return results;
	}

}
