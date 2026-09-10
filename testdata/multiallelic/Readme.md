# snpcheck_multiallelic.bcf

A hand written file for testing `major_minor_alleles_must_be_snps`. It has 6
samples and 5 triallelic sites; at every site two samples carry each of the
three alleles, so all three alleles are observed and can be used as genotypes.

| POS  | REF | ALT   | third allele          | kept |
| ---- | --- | ----- | --------------------- | ---- |
| 1000 | A   | C,G   | SNP                   | yes  |
| 2000 | A   | C,GG  | different length      | no   |
| 3000 | A   | C,ACT | insertion             | no   |
| 4000 | A   | C,\*  | spanning deletion     | no   |
| 5000 | AT  | CT,GT | same length, one base | yes  |

Only the two sites whose alleles are all SNPs may pass the filter. Checking the
two most frequent alleles alone, as was done before, keeps all five.
