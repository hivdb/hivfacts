package edu.stanford.hivdb.hivfacts;

import edu.stanford.hivdb.utilities.ValidationLevel;
import edu.stanford.hivdb.utilities.ValidationMessage;

public enum HIV1ValidationMessage implements ValidationMessage {
	
	NoGeneFound(
		ValidationLevel.CRITICAL,
		"There were no %s genes found, refuse to process."
	),
	FASTAGeneNotAligned(
		ValidationLevel.SEVERE_WARNING,
		"This sequence may also have nucleotides belonging to %s. " +
		"Analysis of this part of the input sequence was suppressed " +
		"due to poor quality, insufficient size or " +
		"improper concatenation of multiple partial sequences."
	),
	FASTAGapTooLong(
		ValidationLevel.CRITICAL,
		"This sequence has critical potentially correctable errors. It has a large sequence gap, " +
		"defined as an insertion or deletion of >30 bps. One possible cause of this error is that " +
		"the input sequence was concatenated from multiple partial sequences. Adding 'N's in place " +
		"of the missing sequence will allow the sequence to be processed."
	),
	FASTASequenceTrimmed(
		ValidationLevel.WARNING,
		"The %s sequence had %d amino acids trimmed from its %s-end due to poor quality."
	),
	FASTASequenceTooShort(
		ValidationLevel.WARNING,
		"The %s sequence contains just %d codons, " +
		"which is not sufficient for a comprehensive interpretation."
	),
	FASTAInvalidNAsRemoved(
		ValidationLevel.NOTE,
		"Non-NA character(s) %s were found and removed from the sequence."
	),
	MultiplePositionsMissing(
		ValidationLevel.WARNING,
		"%d positions were not sequenced or aligned: %s."
	),
	SinglePositionMissing(
		ValidationLevel.NOTE,
		"One position was not sequenced or aligned: %s."
	),
	MultipleDRPsMissing(
		ValidationLevel.WARNING,
		"%d drug-resistance positions " +
		"were not sequenced or aligned: %s."
	),
	SingleDRPMissing(
		ValidationLevel.NOTE,
		"One drug-resistance position " +
		"was not sequenced or aligned: %s."
	),
	NGSMinReadDepthTooLow(
		ValidationLevel.WARNING,
		"You have selected a minimal read-depth of %d. However, " +
		"%d (%.1f%%) position%s in your sequence %s fewer than %d " +
		"reads. Click the ‘Read Coverage’ button to review."
	),
	MultipleStopCodons(
		ValidationLevel.WARNING,
		"There are %d stop codons in %s: %s."
	),
	SingleStopCodon(
		ValidationLevel.WARNING,
		"There is one stop codon in %s: %s."
	),
	MultipleUnusualMutations(
		ValidationLevel.WARNING,
		"There are %d unusual mutations in %s: %s."
	),
	MultipleUnusualMutationsAtDRP(
		ValidationLevel.WARNING,
		"There are %d unusual mutations at drug-resistance positions in %s: %s."
	),
	SingleUnusualMutationAtDRP(
		ValidationLevel.NOTE,
		"There is one unusual mutation at a drug-resistance position in %s: %s."
	),
	SingleAdditionalUnusualMutation(
		ValidationLevel.WARNING,
		"There is one additional unusual mutation in %s: %s."
	),
	MultipleAdditionalUnusualMutations(
		ValidationLevel.SEVERE_WARNING,
		"There are %d additional unusual mutations in %s: %s."
	),
	MultipleApobec(
		ValidationLevel.SEVERE_WARNING,
		"The following %d APOBEC mutations were present in the sequence: %s.%s"
	),
	SingleApobec(
		ValidationLevel.NOTE,
		"This following APOBEC mutation was present in the sequence: %s.%s"
	),
	MultipleApobecAtDRP(
		ValidationLevel.SEVERE_WARNING,
		"There are %d APOBEC-associated mutations at drug-resistance positions: %s."
	),
	SingleApobecAtDRP(
		ValidationLevel.WARNING,
		"There is one APOBEC-associated mutation at a drug-resistance position: %s."
	),
	MultipleUnusualIndelsAndFrameshifts(
		ValidationLevel.SEVERE_WARNING,
		"The %s gene has %d unusual indels and/or frameshifts. " +
		"The indels include %s. The frameshifts include %s."
	),
	MultipleFrameShifts(
		ValidationLevel.SEVERE_WARNING,
		"The %s gene has %d frameshifts: %s."
	),
	MultipleUnusualIndels(
		ValidationLevel.SEVERE_WARNING,
		"The %s gene has %d unusual indels: %s."
	),
	SingleFrameshift(
		ValidationLevel.WARNING,
		"The %s gene has a frameshift: %s."
	),
	SingleUnusualIndel(
		ValidationLevel.WARNING,
		"The %s gene has an unusual indel: %s."
	),
	HIV2(
		ValidationLevel.WARNING,
		"The sequence is from an HIV-2 virus"
	),
	FASTAReverseComplement(
		ValidationLevel.WARNING,
		"This report was derived from the reverse complement of input sequence."
	),
	FASTAUnsequencedRegion(
		ValidationLevel.WARNING,
		"There are %d %s positions located in unsequenced region(s): %s."
	),
	NGSTooManyUnusualMutations(
		ValidationLevel.SEVERE_WARNING,
		"At this threshold (%.1f%%), >= 1.0%% of positions have a " +
		"highly unusual mutation (defined as a prevalence <0.01%% in " +
		"published group M direct PCR sequences). This indicates that " +
		"there may be an unacceptably high risk that some mutations at " +
		"this threshold represent sequence artifacts."
	),
	NGSTooManyApobecMutationsOneApobecDRM(
		ValidationLevel.WARNING,
		"At this threshold (%.1f%%), >=3 positions with signature " +
		"APOBEC mutations. At this threshold, the sequence also contains " +
		"one drug-resistance mutation that could be caused by " +
		"APOBEC-mediated G-to-A hypermutation (%s). This " +
		"DRM therefore should be considered possible sequence artifacts."
	),
	NGSTooManyApobecMutationsMultipleApobecDRMs(
		ValidationLevel.WARNING,
		"At this threshold (%.1f%%), >=3 positions with signature " +
		"APOBEC mutations. At this threshold, the sequence also contains " +
		"%d drug-resistance mutations that could be caused by " +
		"APOBEC-mediated G-to-A hypermutation (%s). These " +
		"DRMs therefore should be considered possible sequence artifacts."
	),
	NGSTooManyApobecMutationsNoApobecDRM(
		ValidationLevel.WARNING,
		"At this threshold (%.1f%%), >=3 positions with signature " +
		"APOBEC mutations. At this threshold, the sequence contains no" +
		"drug-resistance mutations that could be caused by " +
		"APOBEC-mediated G-to-A hypermutation."
	),
	NGSTooLowThreshold(
		ValidationLevel.SEVERE_WARNING,
		"Extensive empirical data, as well as, modeling data suggest that " +
		"the risk of sequence artifact is almost inevitable at thresholds " +
		"below 0.5%% to 1.0%% unless unique molecular identifiers (UMI) are " +
		"used prior to PCR. We display quality control data for the 0.1%% " +
		"and 0.2%% thresholds solely to indicate that choosing a threshold " +
		"that is too low will result in sequence artifacts."
	);
	
	public final ValidationLevel level;
	public final String template;
	
	HIV1ValidationMessage(ValidationLevel level, String template) {
		this.level = level;
		this.template = template;
	}
	
	@Override
	public final ValidationLevel getLevel() {
		return level;
	}
	
	@Override
	public final String getTemplate() {
		return template;
	}
}