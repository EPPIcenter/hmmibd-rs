# source of data file

The original vcf.gz file is from `ftp://ngs.sanger.ac.uk/production/malaria/Resource/34/Pf7_vcf/`

# filtering

BCF file filtering details can be found in the header of the resulting bcf file `pf7_chr1_20samples_maf0.01.bcf`.

# pf7_chr1_20samples_info_ad_dot_ad.bcf

The first 300 records of `pf7_chr1_20samples_maf0.01.bcf`, with two changes used
as a regression test for bcf parsing:

- an `##INFO=<ID=AD,..>` header line is added next to the `##FORMAT=<ID=AD,..>`
  line. Both share one entry of the bcf string dictionary (`IDX=22`), so
  `FORMAT/AD` cannot be found by looking up the dictionary name alone.
- the `FORMAT/AD` value of the first sample is set to `.` in the first three
  records, which must make that sample missing at those sites rather than fail
  the run.
