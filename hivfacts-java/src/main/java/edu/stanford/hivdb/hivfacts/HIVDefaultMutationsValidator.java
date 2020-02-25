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
import java.util.stream.Collectors;

import edu.stanford.hivdb.mutations.Mutation;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationsValidator;
import edu.stanford.hivdb.utilities.ValidationLevel;
import edu.stanford.hivdb.utilities.ValidationResult;

public class HIVDefaultMutationsValidator implements MutationsValidator<HIV> {

	private static final Map<String, ValidationLevel> VALIDATION_RESULT_LEVELS;
	private static final Map<String, String> VALIDATION_RESULT_MESSAGES;

	static {
		Map<String, ValidationLevel> levels = new HashMap<>();
		Map<String, String> messages = new HashMap<>();
		levels.put("severe-warning-too-many-stop-codons", ValidationLevel.SEVERE_WARNING);
		messages.put("severe-warning-too-many-stop-codons", "The submitted mutations contain %d stop codons.");

		levels.put("note-stop-codon", ValidationLevel.WARNING);
		messages.put("note-stop-codon", "The submitted mutations contain %d stop codon.");

		levels.put("much-too-many-unusual-mutations", ValidationLevel.SEVERE_WARNING);
		messages.put("much-too-many-unusual-mutations", "There are %d unusual mutations: %s.");

		levels.put("too-many-unusual-mutations", ValidationLevel.WARNING);
		messages.put("too-many-unusual-mutations", "There are %d unusual mutations: %s.");

		levels.put("some-unusual-mutations", ValidationLevel.NOTE);
		messages.put("some-unusual-mutations", "There are %d unusual mutations: %s.");

		levels.put("unusual-mutation-at-DRP-plural", ValidationLevel.SEVERE_WARNING);
		messages.put("unusual-mutation-at-DRP-plural",
				     "There are %d unusual mutations at drug-resistance positions: %s.");

		levels.put("unusual-mutation-at-DRP", ValidationLevel.WARNING);
		messages.put("unusual-mutation-at-DRP",
				     "There is %d unusual mutation at a drug-resistance position: %s.");

		levels.put("additional-unusual-mutation", ValidationLevel.WARNING);
		messages.put("additional-unusual-mutation", "There is one additional unusual mutation: %s");

		levels.put("additional-unusual-mutations", ValidationLevel.SEVERE_WARNING);
		messages.put("additional-unusual-mutations", "There are %d additional unusual mutations: %s");

		levels.put("severe-APOBEC", ValidationLevel.SEVERE_WARNING);
		messages.put("severe-APOBEC", "The following %d APOBEC muts were present in the sequence.%s");

		levels.put("definite-APOBEC", ValidationLevel.WARNING);
		messages.put("definite-APOBEC", "The following %d APOBEC muts were present in the sequence.%s");

		levels.put("possible-APOBEC-influence", ValidationLevel.NOTE);
		messages.put("possible-APOBEC-influence", "The following %d APOBEC muts were present in the sequence.%s");

		levels.put("multiple-apobec-at-DRP", ValidationLevel.SEVERE_WARNING);
		messages.put("multiple-apobec-at-DRP",
				     "There are %d APOBEC-associated mutations at drug-resistance positions: %s.");

		levels.put("single-apobec-at-DRP", ValidationLevel.WARNING);
		messages.put("single-apobec-at-DRP",
				     "There is %d APOBEC-associated mutation at a drug-resistance position: %s.");

		VALIDATION_RESULT_LEVELS = Collections.unmodifiableMap(levels);
		VALIDATION_RESULT_MESSAGES = Collections.unmodifiableMap(messages);

	}

	@Override
	public List<ValidationResult> validate(MutationSet<HIV> mutations) {
		List<ValidationResult> validationResults = new ArrayList<>();
		validationResults.addAll(validateNoStopCodons(mutations));
		validationResults.addAll(validateNotApobec(mutations));
		validationResults.addAll(validateNoTooManyUnusualMutations(mutations));
		return validationResults;
	}

	private ValidationResult newValidationResult(String key, Object... args) {
		ValidationLevel level = VALIDATION_RESULT_LEVELS.get(key);
		String message = String.format(
			VALIDATION_RESULT_MESSAGES.get(key),
			args);
		return new ValidationResult(level, message);
	}

	private List<ValidationResult> validateNoStopCodons(MutationSet<HIV> mutations) {
		List<ValidationResult> validationResults = new ArrayList<>();
		int numStopCodons = mutations.getStopCodons().size();
		if (numStopCodons > 1) {
			validationResults.add(newValidationResult(
				"severe-warning-too-many-stop-codons",
				numStopCodons));
		} else if (numStopCodons > 0) {
			validationResults.add(newValidationResult(
				"note-stop-codon", numStopCodons));
		}

		return validationResults;
	}


	private List<ValidationResult> validateNoTooManyUnusualMutations(MutationSet<HIV> mutations) {
		List<ValidationResult> validationResults = new ArrayList<>();
		MutationSet<HIV> unusualMutations = mutations.getUnusualMutations();
		int numUnusual = unusualMutations.size();

		MutationSet<HIV> unusualMutAtDRP = unusualMutations.getAtDRPMutations();

		int numUnusualAtDRP = unusualMutAtDRP.size();
		/*if (numUnusual != numUnusualAtDRP && numUnusual > 1) {
			System.out.println("Debug:text:" + text);
			addValidationResult("too-many-unusual-mutations",
					numUnusual, text);
			validated = false;
		} else*/

		if (numUnusualAtDRP > 1) {
			validationResults.add(newValidationResult(
				"unusual-mutation-at-DRP-plural",
				numUnusualAtDRP, unusualMutAtDRP.join(", ", Mutation::getHumanFormatWithGene))
			);
		} else if (numUnusualAtDRP == 1) {
			validationResults.add(newValidationResult(
				"unusual-mutation-at-DRP",
				numUnusualAtDRP,  unusualMutAtDRP.join(", ", Mutation::getHumanFormatWithGene))
			);
		}
		if (numUnusual > 1) {
			if (numUnusualAtDRP == 0) {
				validationResults.add(newValidationResult(
					"too-many-unusual-mutations",
					numUnusual, unusualMutations.join(", ", Mutation::getHumanFormatWithGene))
				);

			} else if (numUnusual - numUnusualAtDRP == 1) {
				MutationSet<HIV> additionalMuts = unusualMutations.subtractsBy(unusualMutAtDRP);
				validationResults.add(newValidationResult(
					"additional-unusual-mutation",
					additionalMuts.join(", ", Mutation::getHumanFormatWithGene))
				);

			} else if (numUnusual -numUnusualAtDRP > 1) {
				int numAdditionalUnusual = numUnusual - numUnusualAtDRP;
				MutationSet<HIV> additionalMuts = unusualMutations.subtractsBy(unusualMutAtDRP);
				validationResults.add(newValidationResult(
					"additional-unusual-mutations", numAdditionalUnusual,
					additionalMuts.join(", ", Mutation::getHumanFormatWithGene))
				);
			} else {
				// numUnusual == numUnusualAtDRP && numUnusual > 1
				// No  warning  necessary as it will be included in the unusual-mutation-at-DRP-plural warning
			}
		}

		return validationResults;

	}

	private List<ValidationResult> validateNotApobec(MutationSet<HIV> mutations) {
		MutationSet<HIV> apobecs = mutations.getApobecMutations();
		MutationSet<HIV> apobecDRMs = mutations.getApobecDRMs();
		List<ValidationResult> results = new ArrayList<>();
		int numApobecMuts = apobecs.size();
		int numApobecDRMs = apobecDRMs.size();
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
				)).collect(Collectors.joining(","))
			);
		}

		if (numApobecMuts > 3) {
			results.add(newValidationResult("severe-APOBEC", numApobecMuts, extraCmt));
		} else if (numApobecMuts > 1) {
			results.add(newValidationResult("definite-APOBEC", numApobecMuts, extraCmt));
		} else if (numApobecMuts == 1) {
			results.add(newValidationResult("possible-APOBEC-influence", numApobecMuts, extraCmt));
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

}