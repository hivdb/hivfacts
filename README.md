# HIV-AAPcnt

[![Travis CI](https://api.travis-ci.org/hivdb/hiv-aapcnt.svg?branch=master)](https://travis-ci.org/hivdb/hiv-aapcnt)
[![codecov](https://codecov.io/gh/hivdb/hiv-aapcnt/branch/master/graph/badge.svg)](https://codecov.io/gh/hivdb/hiv-aapcnt)

Amino acid prevalence data of HIV-1 pol.

## Steps and Criteria for Classifying Unusual Mutations

First we calculate the prevalence of an amino acid at a position:

1. Same amino acid at the same position from a single person is only
   counted once. (`PatientPerMutCount`)
2. Multiple Amino acids at the same position from a single sample are not
   counted. (`excludeMixtures`)
3. Stop codons are not counted. (`excludeStopCodons`)

In general, an amino acid at a position with its prevalence < 0.01% is
considered unusual, with the following exceptions:

1. It is possessed with drug resistance score in HIVDB algorithm;
2. It is considered as a signature APOBEC-mediated hypermutation.
