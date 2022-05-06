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
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationsValidator;
import edu.stanford.hivdb.utilities.MyStringUtils;
import edu.stanford.hivdb.utilities.ValidationLevel;
import edu.stanford.hivdb.utilities.ValidationResult;
import edu.stanford.hivdb.viruses.Gene;

public class HIVDefaultMutationsValidator implements MutationsValidator<HIV> {

	@Override
	public List<ValidationResult> validate(MutationSet<HIV> mutations, Collection<String> includeGenes) {
		List<ValidationResult> validationResults = new ArrayList<>();
		validationResults.addAll(validateNoStopCodons(mutations, includeGenes));
		validationResults.addAll(validateNotApobec(mutations, includeGenes));
		validationResults.addAll(validateNoTooManyUnusualMutations(mutations, includeGenes));
		return validationResults;
	}

	private static List<ValidationResult> validateNoStopCodons(
		MutationSet<HIV> mutations,
		Collection<String> includeGenes
	) {
		List<ValidationResult> validationResults = new ArrayList<>();
		MutationSet<HIV> stopCodons = mutations
			.getStopCodons()
			.filterBy(mut -> includeGenes.contains(mut.getAbstractGene()));
		for (Map.Entry<Gene<HIV>, MutationSet<HIV>> entry : stopCodons.groupByGene().entrySet()) {
			String geneText = entry.getKey().getAbstractGene();
			MutationSet<HIV> geneStopCodons = entry.getValue();
			int numGeneStopCodons = geneStopCodons.size();
			String geneStopText = geneStopCodons.join(", ", Mutation::getHumanFormatWithAbstractGene);
			if (numGeneStopCodons > 1) {
				validationResults.add(HIV1ValidationMessage.MultipleStopCodons.formatWithLevel(
					ValidationLevel.SEVERE_WARNING,
					numGeneStopCodons,
					geneText,
					geneStopText
				));
			} else if (numGeneStopCodons > 0) {
				validationResults.add(HIV1ValidationMessage.SingleStopCodon.formatWithLevel(
					ValidationLevel.WARNING,
					geneText,
					geneStopText
				));
			}
		}
		
		return validationResults;
	}


	private static List<ValidationResult> validateNoTooManyUnusualMutations(
		MutationSet<HIV> mutations,
		Collection<String> includeGenes
	) {
		List<ValidationResult> validationResults = new ArrayList<>();
		MutationSet<HIV> unusualMuts = mutations
			.getUnusualMutations()
			.filterBy(mut -> includeGenes.contains(mut.getAbstractGene()));

		for (Map.Entry<Gene<HIV>, MutationSet<HIV>> entry : unusualMuts.groupByGene().entrySet()) {
			String geneText = entry.getKey().getAbstractGene();
			MutationSet<HIV> geneUnusualMuts = entry.getValue();
			int numGeneUnusual = geneUnusualMuts.size();

			MutationSet<HIV> geneUnusualMutsAtDRP = geneUnusualMuts.getAtDRPMutations();
			int numGeneUnusualAtDRP = geneUnusualMutsAtDRP.size();

			if (numGeneUnusualAtDRP > 1) {
				validationResults.add(HIV1ValidationMessage.MultipleUnusualMutationsAtDRP.formatWithLevel(
					ValidationLevel.SEVERE_WARNING,
					numGeneUnusualAtDRP,
					geneText,
					geneUnusualMutsAtDRP.join(", "))
				);
			} else if (numGeneUnusualAtDRP == 1) {
				validationResults.add(HIV1ValidationMessage.SingleUnusualMutationAtDRP.formatWithLevel(
						ValidationLevel.WARNING,
						geneText,
						geneUnusualMutsAtDRP.join(", "))
				);
			}
			if (numGeneUnusual > 1) {
				String geneMutText = geneUnusualMuts.join(", ");
				if (numGeneUnusualAtDRP == 0) {
					validationResults.add(HIV1ValidationMessage.MultipleUnusualMutations.formatWithLevel(
						ValidationLevel.WARNING,
						numGeneUnusual,
						geneText,
						geneMutText
					));

				} else if (numGeneUnusual - numGeneUnusualAtDRP == 1) {
					MutationSet<HIV> additionalMuts = geneUnusualMuts.subtractsBy(geneUnusualMutsAtDRP);
					validationResults.add(HIV1ValidationMessage.SingleAdditionalUnusualMutation.formatWithLevel(
						ValidationLevel.WARNING,
						geneText,
						additionalMuts.join(", ")
					));
					

				} else if (numGeneUnusual - numGeneUnusualAtDRP > 1) {
					int numAdditionalUnusual = numGeneUnusual - numGeneUnusualAtDRP;
					MutationSet<HIV> additionalMuts = geneUnusualMuts.subtractsBy(geneUnusualMutsAtDRP);
					validationResults.add(HIV1ValidationMessage.MultipleAdditionalUnusualMutations.formatWithLevel(
						ValidationLevel.SEVERE_WARNING,
						numAdditionalUnusual,
						geneText,
						additionalMuts.join(", ")
					));
				} else {
					// numGeneUnusual == numGeneUnusualAtDRP && numGeneUnusual > 1
					// No warning necessary as it will be included in the MultipleUnusualMutationsAtDRP warning
				}
			}
		}

		return validationResults;

	}

	private static List<ValidationResult> validateNotApobec(
		MutationSet<HIV> mutations,
		Collection<String> includeGenes
	) {
		MutationSet<HIV> apobecs = mutations
			.getApobecMutations()
			.filterBy(mut -> includeGenes.contains(mut.getAbstractGene()));
		MutationSet<HIV> apobecDRMs = mutations
			.getApobecDRMs()
			.filterBy(mut -> includeGenes.contains(mut.getAbstractGene()));
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
					MyStringUtils.andListFormat(e.getValue())
				)).collect(Collectors.joining("; "))
			);
		}

		if (numApobecMuts > 3) {
			results.add(HIV1ValidationMessage.MultipleApobec.formatWithLevel(
				ValidationLevel.SEVERE_WARNING,
				numApobecMuts,
				apobecMutsText,
				extraCmt
			));
		} else if (numApobecMuts > 1) {
			results.add(HIV1ValidationMessage.MultipleApobec.formatWithLevel(
				ValidationLevel.WARNING,
				numApobecMuts,
				apobecMutsText,
				extraCmt
			));
		} else if (numApobecMuts == 1) {
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

}