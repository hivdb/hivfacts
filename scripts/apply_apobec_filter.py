#! /usr/bin/env python3

import sys
import csv
import json

THRESHOLD_TARGET_AA_PCNT = 0.5
THRESHOLD_TARGET_CODON_PCNT = 0.01
THRESHOLD_TARGET_CODON_COUNT = 1000

THRESHOLD_MUTATED_CODON_PCNT_MAJOR_SUBTYPE = 0.01
THRESHOLD_MUTATED_CODON_PCNT_OTHER_SUBTYPE = 0.02

THRESHOLD_MAX_MUTATED_AA_PCNT = 0.005
THRESHOLD_STOP_RATE_FOR_ONE_IN_THOUSAND = 0.4
THRESHOLD_STOP_RATE_FOR_ONE_IN_TEN_THOUSANDS = 0.3


def main():
    if len(sys.argv) != 7:
        print('Usage: {} <SPECIES> <ALL_HM_JSON> <APOBEC_JSON> '
              '<APOBEC_DRM_JSON> <APOBEC_CSV> <APOBEC_DRM_CSV>',
              file=sys.stderr)
        exit(1)
    species, hm_input, hm_out, hmdrm_out, hm_csv, hmdrm_csv = sys.argv[1:]
    with open(hm_input) as fp:
        hmdata = json.load(fp)

    step_1a = [r for r in hmdata
               if r['target_aa_prevalence'] > THRESHOLD_TARGET_AA_PCNT]
    step_1b = [r for r in step_1a
               if (r['target_codon_prevalence'] or 0) >
               THRESHOLD_TARGET_CODON_PCNT
               or
               (r['target_codon_prevalence'] or 0) *
               r['pos_total'] > THRESHOLD_TARGET_CODON_COUNT]
    print('Removed by step 1b:', len(step_1b))
    step_3a0 = [r for r in step_1b if not r['is_drm']]
    print('Removed by step 3a0:', len(step_1b) - len(step_3a0))
    step_3a1 = [r for r in step_3a0
                if r['mutated_aa_prevalence'] <=
                THRESHOLD_MAX_MUTATED_AA_PCNT]
    print('Removed by step 3a1:', len(step_3a0) - len(step_3a1))
    if species == 'HIV1':
        step_3a2 = [r for r in step_3a1
                    if max(
                        r['mutated_codon_prevalence_A'] or 0,
                        r['mutated_codon_prevalence_B'] or 0,
                        r['mutated_codon_prevalence_C'] or 0,
                        r['mutated_codon_prevalence_CRF01_AE'] or 0,
                        r['mutated_codon_prevalence_CRF02_AG'] or 0,
                    ) <=
                    THRESHOLD_MUTATED_CODON_PCNT_MAJOR_SUBTYPE]
    else:
        step_3a2 = [r for r in step_3a1
                    if max(
                        r['mutated_codon_prevalence_Group A'] or 0,
                        r['mutated_codon_prevalence_Group B'] or 0,
                    ) <=
                    THRESHOLD_MUTATED_CODON_PCNT_MAJOR_SUBTYPE]
    print('Removed by step 3a2:', len(step_3a1) - len(step_3a2))
    step_3a3 = [r for r in step_3a2
                if (r['mutated_codon_prevalence_others'] or 0) <=
                THRESHOLD_MUTATED_CODON_PCNT_OTHER_SUBTYPE]
    print('Removed by step 3a3:', len(step_3a2) - len(step_3a3))
    step_3b1 = [r for r in step_3a3
                if r['mutated_aa_prevalence'] <= 0.001 or
                r['pcnt_stops_all'] > THRESHOLD_STOP_RATE_FOR_ONE_IN_THOUSAND]
    print('Removed by step 3b1:', len(step_3a3) - len(step_3b1))
    step_3b2 = [r for r in step_3b1
                if r['mutated_aa_prevalence'] <= 0.0001 or
                r['pcnt_stops_all'] >
                THRESHOLD_STOP_RATE_FOR_ONE_IN_TEN_THOUSANDS]
    print('Removed by step 3b2:', len(step_3b1) - len(step_3b2))
    print('APOBECs:', len(step_3b2))
    apobec_drms = [r for r in step_1a if r['is_drm']]
    print('APOBEC-context DRMs:', len(apobec_drms))

    with open(hm_out, 'w') as fp:
        json.dump([{
            'gene': r['gene'],
            'position': r['codon_pos'],
            'aa': r['mutated_aa']
        } for r in step_3b2], fp, indent=2)

    with open(hmdrm_out, 'w') as fp:
        json.dump([{
            'gene': r['gene'],
            'position': r['codon_pos'],
            'aa': r['mutated_aa']
        } for r in apobec_drms], fp, indent=2)

    header = list(step_3b2[0].keys())
    with open(hm_csv, 'w') as fp:
        writer = csv.DictWriter(fp, header)
        writer.writeheader()
        for row in step_3b2:
            row['target_aas'] = ''.join(row['target_aas'])
            row['codon_pairs'] = ', '.join(
                '{0}{2} => {1}{2} ({3})'.format(l, r, c, t)
                for l, r, c, t in row['codon_pairs'])
        writer.writerows(step_3b2)

    with open(hmdrm_csv, 'w') as fp:
        writer = csv.DictWriter(fp, header)
        writer.writeheader()
        for row in apobec_drms:
            row['target_aas'] = ''.join(row['target_aas'])
            row['codon_pairs'] = ', '.join(
                '{0}{2} => {1}{2} ({3})'.format(l, r, c, t)
                for l, r, c, t in row['codon_pairs'])
        writer.writerows(apobec_drms)


if __name__ == '__main__':
    main()
