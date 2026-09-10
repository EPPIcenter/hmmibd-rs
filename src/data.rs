use itertools::{EitherOrBoth, Itertools};
use smallvec::SmallVec;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read},
    path::Path,
    sync::{Arc, Mutex},
};

use crate::{
    args::Arguments,
    bcf::{self, BcfFilterArgs, BcfGenotype},
    genome::Genome,
    matrix::*,
    samples::{self, Samples},
    sites::{self, SiteInfoRaw, Sites},
};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("error in processing bcf file: {0:#?}")]
    Bcf(#[from] crate::bcf::Error),
    #[error("io error, source {source:?}, file: {file:?}")]
    Io {
        source: std::io::Error,
        file: Option<String>,
    },
    // #[error("custom: {0:?}")]
    // Custom(String),
    #[error("{0:?}")]
    ParseLineError(#[from] ParseLineError),
    #[error("irow={irow}, iallele={iallele}, value={value:?}")]
    MissingFrequency {
        irow: usize,
        iallele: usize,
        value: Option<f64>,
    },
    #[error("EmptyIterator")]
    EmptyIterator,

    #[error("lockerror: {0}")]
    LockError(&'static str),

    #[error("sites error: {0:?}")]
    Site(#[from] sites::Error),

    #[error("sample error: {0:?}")]
    Sample(#[from] samples::Error),

    #[error("sample error: {0:?}")]
    PartialCmpIsNone(&'static str),

    #[error("toml serialize error: {0:?}")]
    TomlSerializeError(#[from] toml::ser::Error),

    #[error("cli argument error: {0:?}")]
    CliArgError(&'static str),

    #[error("cli argument error: output prefix cannot be inferred")]
    OutputPrefixCantBeInferred,

    #[error("frequency inference error: genotype missing across all samples at a given site")]
    FreqInferenceZeroNonMissGenotype,

    #[error(
        "frequency inference error: genotype missing across all samples listed in \
        {file} at site {chrname}:{pos}"
    )]
    FreqInferenceZeroNonMissGenotypeInSubset {
        file: String,
        chrname: String,
        pos: u32,
    },

    #[error("none of the sample ids listed in {0} is among the analyzed samples")]
    FreqSamplesNotFound(String),

    #[error(
        "sample {sample} listed in {file} does not belong to the population the \
        option applies to"
    )]
    FreqSampleWrongPopulation { sample: String, file: String },

    #[error(
        "the two populations share no site; check that --data-file1 and \
        --data-file2 are aligned to the same reference and that the bcf \
        filtering options are not too stringent"
    )]
    NoSharedSitesBetweenPopulations,
}

#[derive(Debug, thiserror::Error)]
pub enum ParseLineError {
    #[error("Cannot read column {0}")]
    ReadColumnError(&'static str),
    #[error("Cannot parse column {0}")]
    ParseColumnError(&'static str),
    #[error("Cannot read line")]
    ReadLineError,
}

pub struct InputData {
    /// arguments
    pub args: Arguments,
    /// genotype matrix for population 1
    pub geno: Matrix<u8>,
    /// an array of number of unique alleles, length number sites
    pub nall: Vec<u8>,
    pub majall: Vec<u8>,
    /// genotype matrix for population 2
    pub freq1: Matrix<f64>,
    /// allele frequency matrix for population 2
    pub freq2: Option<Matrix<f64>>,
    pub sites: Sites,
    pub genome: Genome,
    pub samples: Samples,
    pub pairs: Vec<(u32, u32)>,
}

impl InputData {
    /// The first element of return value is always a chunk index in the samples in pop1
    /// The 2nd element of return value is a chunk index that can be in pop1 or pop2.
    /// the chunk index for pop2 is always no less than max chunk index of pop1 which can be used
    /// to determine which population the chunk is from
    pub fn get_chunk_pairs(&self) -> Vec<(u32, u32)> {
        let mut chunk_pairs = vec![];
        let mut num_chunk1 = self.samples.pop1_nsam() / self.args.par_chunk_size;
        if self.samples.pop1_nsam() % self.args.par_chunk_size > 0 {
            num_chunk1 += 1;
        }
        if self.freq2.is_some() {
            // two different populations
            assert!(self.samples.pop2_nsam() > 0);
            let mut num_chunk2 = self.samples.pop2_nsam() / self.args.par_chunk_size;
            if self.samples.pop2_nsam() % self.args.par_chunk_size > 0 {
                num_chunk2 += 1;
            }
            for i in 0..num_chunk1 {
                for j in 0..num_chunk2 {
                    let j = num_chunk1 + j;
                    chunk_pairs.push((i, j));
                }
            }
        } else {
            // same population
            for i in 0..num_chunk1 {
                for j in i..num_chunk1 {
                    chunk_pairs.push((i, j));
                }
            }
        }
        chunk_pairs
    }

    pub fn clone_inputdata_for_chunkpair(&self, chunkpair: (u32, u32)) -> Result<Self, Error> {
        let (ichunk, jchunk) = chunkpair;

        // num of chunks for samples.s1()
        let num_chunk1 = {
            let a = self.samples.pop1_nsam();
            let b = self.args.par_chunk_size;
            let mut c = a / b;
            let d = a % b;
            if d > 0 {
                c += 1;
            }
            c
        };

        // chunk1 start,end rows/samples, len
        let (ichunk_start, ichunk_end, ichunk_len) = {
            let ichunk_start = ichunk * self.args.par_chunk_size;
            let mut ichunk_end = (1 + ichunk) * self.args.par_chunk_size;
            if ichunk_end > self.samples.pop1_nsam() {
                ichunk_end = self.samples.pop1_nsam();
            }
            let ichunk_len = ichunk_end - ichunk_start;
            (ichunk_start, ichunk_end, ichunk_len)
        };

        // chunk2  start, end rows/samples, len
        let (jchunk_start, jchunk_end, jchunk_len) = {
            if jchunk < num_chunk1 {
                let jchunk_start = jchunk * self.args.par_chunk_size;
                let mut jchunk_end = (1 + jchunk) * self.args.par_chunk_size;
                if jchunk_end > self.samples.pop1_nsam() {
                    jchunk_end = self.samples.pop1_nsam();
                }
                let jchunk_len = jchunk_end - jchunk_start;
                (jchunk_start, jchunk_end, jchunk_len)
            } else {
                assert!(self.samples.pop2_nsam() > 0);
                let jjchunk = jchunk - num_chunk1;
                let jchunk_start = self.samples.pop1_nsam() + jjchunk * self.args.par_chunk_size;
                let mut jchunk_end =
                    self.samples.pop1_nsam() + (1 + jjchunk) * self.args.par_chunk_size;
                if jchunk_end > self.samples.pop1_nsam() + self.samples.pop2_nsam() {
                    jchunk_end = self.samples.pop1_nsam() + self.samples.pop2_nsam();
                }
                (jchunk_start, jchunk_end, jchunk_end - jchunk_start)
            }
        };
        let nsites = self.geno.get_ncols();

        let args = self.args.clone();
        let nall = self.nall.clone();
        let majall = self.majall.clone();
        let freq1 = self.freq1.clone();
        let geno = {
            let mut nrows = ichunk_len as usize;
            let ncols = nsites;
            let mut s = ichunk_start as usize * nsites;
            let mut e = ichunk_end as usize * nsites;
            let mut data = self.geno.as_slice()[s..e].to_vec();
            if ichunk != jchunk {
                s = jchunk_start as usize * nsites;
                e = jchunk_end as usize * nsites;
                data.extend_from_slice(&self.geno.as_slice()[s..e]);
                nrows += jchunk_len as usize;
            }
            Matrix::<u8>::from_shape_vec(nrows, ncols, data)
        };
        let freq2 = match (ichunk == jchunk, jchunk < num_chunk1) {
            (true, _) => None,
            (false, true) => Some(self.freq1.clone()),
            (false, false) => self.freq2.clone(),
        };
        let samples = {
            let pop1_slice = &self.samples.v()[(ichunk_start as usize)..(ichunk_end as usize)];
            if ichunk == jchunk {
                Samples::from_slice(pop1_slice, &[])
            } else {
                let pop2_slice = &self.samples.v()[(jchunk_start as usize)..(jchunk_end as usize)];
                Samples::from_slice(pop1_slice, pop2_slice)
            }
        };
        let sites = self.sites.clone();
        let genome = self.genome.clone();
        let pairs = Self::get_valid_pair_file(&args, &samples)?;
        Ok(Self {
            args,
            geno,
            nall,
            majall,
            freq1,
            freq2,
            sites,
            genome,
            samples,
            pairs,
        })
    }

    pub fn from_args(args: &Arguments) -> Result<Self, Error> {
        let bcf_gt = if args.from_bin {
            Some(BcfGenotype::load_from_file(&args.data_file1)?)
        } else if args.from_bcf {
            let bcf_filter_args = match args.bcf_filter_config.as_ref() {
                Some(dom_gt_config_path) => BcfFilterArgs::new_from_toml_file(dom_gt_config_path),
                None => {
                    let config = BcfFilterArgs::new_from_builtin()?;
                    std::fs::write("tmp_bcf_filter_config.toml", toml::to_string(&config)?)
                        .map_err(|e| Error::Io {
                            source: e,
                            file: args.bcf_filter_config.to_owned(),
                        })?;
                    eprintln!(concat!(
                        "WARN: bcf_filter_config not specified, a builtin configuration is used",
                        " and is written to 'tmp_bcf_filter_config.toml'",
                    ));
                    Ok(config)
                }
            }?;
            let bcf_gt = BcfGenotype::new_from_processing_bcf(
                &args.bcf_read_mode,
                &bcf_filter_args,
                &args.data_file1,
            )?;

            // check if --bcf-to-bin-file-by-chromosome
            if args.bcf_to_bin_file_by_chromosome {
                bcf_gt.split_chromosomes_into_files(
                    args.output
                        .as_ref()
                        .ok_or(Error::OutputPrefixCantBeInferred)?,
                )?;
            }
            // check if --bcf-to-bin-file
            if args.bcf_to_bin_file {
                let output_prefix = args
                    .output
                    .as_ref()
                    .ok_or(Error::OutputPrefixCantBeInferred)?;
                if let Some(parent) = std::path::Path::new(output_prefix).parent() {
                    std::fs::create_dir_all(parent).map_err(|e| Error::Io {
                        source: e,
                        file: Some(format!("{parent:?}")),
                    })?;
                }
                let output = format!("{}.bin", output_prefix);
                bcf_gt.save_to_file(&output)?;
            }
            if args.bcf_to_bin_file_by_chromosome || args.bcf_to_bin_file {
                // write sample name list
                let output_prefix = args
                    .output
                    .as_ref()
                    .ok_or(Error::OutputPrefixCantBeInferred)?;
                let output = format!("{}.samples", output_prefix);

                std::fs::write(&output, bcf_gt.get_samples().join("\n")).map_err(|e| {
                    Error::Io {
                        source: e,
                        file: Some(output),
                    }
                })?;

                eprintln!(
                    "WARN: --bcf-to-bin-file-xxx is/are specified; \
                    bin file(s) have been written; hmm inference is skipped"
                );
                std::process::exit(0);
            }

            Some(bcf_gt)
        } else {
            None
        };

        // genotype of a second population, when `-I/--data-file2` is used
        // together with `--from-bcf` or `--from-bin`
        let bcf_gt2 = match (bcf_gt.as_ref(), args.data_file2.as_ref()) {
            (Some(_), Some(data_file2)) if args.from_bin => {
                Some(BcfGenotype::load_from_file(data_file2)?)
            }
            (Some(_), Some(data_file2)) => {
                let bcf_filter_args = match args.bcf_filter_config.as_ref() {
                    Some(dom_gt_config_path) => {
                        BcfFilterArgs::new_from_toml_file(dom_gt_config_path)?
                    }
                    None => BcfFilterArgs::new_from_builtin()?,
                };
                Some(BcfGenotype::new_from_processing_bcf(
                    &args.bcf_read_mode,
                    &bcf_filter_args,
                    data_file2,
                )?)
            }
            _ => None,
        };

        let valid_samples = Samples::from_args(args, bcf_gt.as_ref(), bcf_gt2.as_ref())?;
        let min_snp_sep = args.min_snp_sep;

        let (geno1, geno2, sitesinfo) = match args.from_bcf || args.from_bin {
            true => {
                let dgt = bcf_gt.ok_or(bcf::Error::GenotypeEmpty {
                    file: file!(),
                    line: line!(),
                })?;
                match bcf_gt2 {
                    None => {
                        let (geno1, sitesinfo) =
                            Self::read_data_dominant_genotype(dgt, &valid_samples, min_snp_sep)?;
                        (geno1, None, sitesinfo)
                    }
                    Some(dgt2) => {
                        // each bcf file is filtered on its own, so the two site
                        // sets generally differ; the sites are thinned by
                        // `--min-snp-sep` only after intersecting them so that
                        // the kept sites do not depend on which sites a single
                        // population happens to have
                        let (geno1, sitesinfo1) =
                            Self::read_data_dominant_genotype(dgt, &valid_samples, 0)?;
                        let (geno2, sitesinfo2) =
                            Self::read_data_dominant_genotype(dgt2, &valid_samples, 0)?;
                        let (geno1, geno2, sitesinfo) = Self::intersect_sites_of_two_populations(
                            geno1,
                            sitesinfo1,
                            geno2,
                            sitesinfo2,
                            min_snp_sep,
                        )?;
                        (geno1, Some(geno2), sitesinfo)
                    }
                }
            }
            false => {
                let (geno1, sitesinfo) =
                    Self::read_data_hmmibd_format(&args.data_file1, &valid_samples, min_snp_sep)?;
                let geno2 = match args.data_file2.as_ref() {
                    Some(data_file2) => {
                        let (geno2, sites2) =
                            Self::read_data_hmmibd_format(data_file2, &valid_samples, min_snp_sep)?;
                        assert_eq!(&sitesinfo, &sites2);
                        Some(geno2)
                    }
                    None => None,
                };
                (geno1, geno2, sitesinfo)
            }
        };

        // samples the allele frequencies are calculated from, when restricted
        // by `--freq-samples1`/`--freq-samples2`
        let freq_cols1 = match (args.freq_file1.as_ref(), args.freq_samples1.as_ref()) {
            (None, Some(freq_samples1)) => Some(Self::get_freq_sample_cols(
                freq_samples1,
                &valid_samples,
                0,
                geno1.get_ncols(),
            )?),
            _ => None,
        };
        let freq_cols2 = match (
            geno2.as_ref(),
            args.freq_file2.as_ref(),
            args.freq_samples2.as_ref(),
        ) {
            (Some(geno2), None, Some(freq_samples2)) => Some(Self::get_freq_sample_cols(
                freq_samples2,
                &valid_samples,
                valid_samples.pop1_nsam() as usize,
                geno2.get_ncols(),
            )?),
            _ => None,
        };
        let (mut geno1, geno2, sitesinfo) = Self::drop_sites_without_freq_samples(
            geno1,
            geno2,
            sitesinfo,
            freq_cols1.as_deref(),
            freq_cols2.as_deref(),
        );

        // create the genome object and site object
        let (sites, genome) = sitesinfo.into_sites_and_genome(&args.rec_args)?;

        let nsam1_valid = geno1.get_ncols();
        let nsam2_valid = geno2.as_ref().map(|geno2| geno2.get_ncols());
        let freq1 = match (args.freq_file1.as_ref(), freq_cols1.as_ref()) {
            (Some(freq_file1), _) => Self::read_freq_file(freq_file1, args, &genome, &sites)?,
            (None, Some(cols)) => Self::infer_freq_from_data_subset(
                &geno1,
                &sites,
                &genome,
                args,
                cols,
                args.freq_samples1.as_deref().unwrap_or(""),
            )?,
            (None, None) => Self::infer_freq_from_data(&geno1, &sites, args)?,
        };

        let freq2 = match geno2.as_ref() {
            None => None,
            Some(geno2) => match (args.freq_file2.as_ref(), freq_cols2.as_ref()) {
                (Some(freq_file2), _) => {
                    Some(Self::read_freq_file(freq_file2, args, &genome, &sites)?)
                }
                (None, Some(cols)) => Some(Self::infer_freq_from_data_subset(
                    geno2,
                    &sites,
                    &genome,
                    args,
                    cols,
                    args.freq_samples2.as_deref().unwrap_or(""),
                )?),
                (None, None) => Some(Self::infer_freq_from_data(geno2, &sites, args)?),
            },
        };
        let pairs = Self::get_valid_pair_file(args, &valid_samples)?;
        let nall = Self::get_nall(&freq1, freq2.as_ref())?;
        let majall = Self::get_major_all(&freq1, freq2.as_ref(), nsam1_valid, nsam2_valid)?;

        // sample oriented
        geno1.transpose();
        if let Some(mut geno2) = geno2 {
            geno2.transpose();
            geno1.merge(&geno2);
        }

        Ok(Self {
            args: args.clone(),
            geno: geno1,
            nall,
            majall,
            freq1,
            freq2,
            samples: valid_samples,
            sites,
            genome,
            pairs,
        })
    }
    // fn read_data_file_with_ginfo_and_gmap(
    //     data_file: impl AsRef<Path>,
    //     valid_samples: &Samples,
    //     min_snp_sep: u32,
    //     ginfo: &GenomeInfo,
    //     gmap: &GeneticMap,
    //     genome: &Genome,
    // ) -> (Matrix<u8>, Sites) {
    //     let mut sites = Sites::new();

    //     let mut line = String::with_capacity(100000);
    //     let mut f = std::fs::File::open(data_file.as_ref())
    //         .map(BufReader::new)
    //         .unwrap();

    //     // get all sample names
    //     f.read_line(&mut line).unwrap();
    //     let samples: Vec<_> = line
    //         .trim()
    //         .split('\t')
    //         .skip(2)
    //         .map(|x| x.to_owned())
    //         .collect();
    //     line.clear();

    //     let n_valid_samples = samples
    //         .iter()
    //         .filter(|s| valid_samples.m().contains_key(*s))
    //         .count();

    //     let mut geno1 = MatrixBuilder::<u8>::new(n_valid_samples);
    //     let mut last_chrname = String::new();
    //     let mut last_pos = 0;
    //     while f.read_line(&mut line).unwrap() != 0 {
    //         let mut fields = line.trim().split("\t");
    //         let chrname = fields.next().unwrap();
    //         let pos: u32 = fields.next().unwrap().parse().unwrap();

    //         // println!("pos: {pos}, last_pos: {last_pos}");
    //         if (chrname == last_chrname) && (last_pos + min_snp_sep > pos) {
    //             line.clear();
    //             continue;
    //         } else {
    //             last_chrname.clear();
    //             last_chrname.push_str(chrname);
    //             last_pos = pos;
    //         }

    //         let chrid = ginfo.idx[chrname];
    //         let gw_pos = ginfo.to_gw_pos(chrid, pos);
    //         let gw_pos_cm = gmap.get_cm(gw_pos);

    //         sites.add(gw_pos, gw_pos_cm);

    //         let mut cnt = 0;
    //         fields.enumerate().for_each(|(i, field)| {
    //             if valid_samples.m().contains_key(&samples[i]) {
    //                 let allel = match field.parse::<i8>().unwrap() {
    //                     -1 => None,
    //                     x => Some(x as u8),
    //                 };
    //                 geno1.push(allel);
    //                 cnt += 1;
    //             }
    //         });
    //         assert_eq!(cnt, n_valid_samples);
    //         line.clear();
    //     }

    //     let geno1 = geno1.finish();
    //     sites.finish(&genome);

    //     (geno1, sites)
    // }

    fn read_data_dominant_genotype(
        dgt: BcfGenotype,
        valid_samples: &Samples,
        min_snp_sep: u32,
    ) -> Result<(Matrix<u8>, SiteInfoRaw), Error> {
        // check seleted_samples are correct
        dgt.into_genotype_siteinfo(valid_samples, min_snp_sep)
            .map_err(|e| e.into())
    }

    fn read_data_hmmibd_format(
        data_file: impl AsRef<Path>,
        valid_samples: &Samples,
        min_snp_sep: u32,
    ) -> Result<(Matrix<u8>, SiteInfoRaw), Error> {
        let mut siteinfo = SiteInfoRaw::new();

        let mut line = String::with_capacity(100000);
        let mut f = std::fs::File::open(data_file.as_ref())
            .map(BufReader::new)
            .map_err(|source| Error::Io {
                source,
                file: Some(data_file.as_ref().to_string_lossy().to_string()),
            })?;

        // get all sample names
        f.read_line(&mut line).map_err(|e| Error::Io {
            source: e,
            file: None,
        })?;
        let samples: Vec<_> = line
            .trim()
            .split('\t')
            .skip(2)
            .map(|x| x.to_owned())
            .collect();
        line.clear();

        let valid_sample_col: Vec<usize> = samples
            .iter()
            .enumerate()
            .filter(|(_i, s)| valid_samples.m().contains_key(*s))
            .map(|(i, _)| i)
            .collect();
        let n_valid_samples = valid_sample_col.len();

        let mut geno = MatrixBuilder::<u8>::new(n_valid_samples);
        let mut last_chrname = String::new();
        let mut last_pos = 0;

        while f
            .read_line(&mut line)
            .map_err(|source| Error::Io { source, file: None })?
            != 0
        {
            let mut fields = line.trim().split("\t");
            let chrname = fields
                .next()
                .ok_or(ParseLineError::ReadColumnError("chrname"))?;
            let pos: u32 = fields
                .next()
                .ok_or(ParseLineError::ReadColumnError("pos"))?
                .parse()
                .map_err(|_| ParseLineError::ParseColumnError("pos"))?;

            // println!("pos: {pos}, last_pos: {last_pos}");
            if chrname == last_chrname {
                if last_pos + min_snp_sep > pos {
                    line.clear();
                    continue;
                } else {
                    last_chrname.clear();
                    last_chrname.push_str(chrname);
                    last_pos = pos;
                }
            } else {
                last_chrname.clear();
                last_chrname.push_str(chrname);
                last_pos = pos;
                siteinfo.add_chr_name(chrname);
            }

            siteinfo.add_chr_idx(chrname);
            siteinfo.add_chr_pos(pos);

            let mut cnt = 0;

            // use merge_join_by to only parse columns that contain selected samples
            fields
                .enumerate()
                .merge_join_by(valid_sample_col.iter(), |a, b| a.0.cmp(b))
                .try_for_each(|mergeby_res| -> Result<(), Error> {
                    if let EitherOrBoth::Both((_, field), _) = mergeby_res {
                        let allel = match field
                            .parse::<i8>()
                            .map_err(|_| ParseLineError::ParseColumnError("alleles"))?
                        {
                            -1 => None,
                            x => Some(x as u8),
                        };
                        geno.push(allel);
                        cnt += 1;
                    }
                    Ok(())
                })?;
            assert_eq!(cnt, n_valid_samples);
            line.clear();
        }

        let geno = geno.finish();
        Ok((geno, siteinfo))
    }

    // fn read_data_file(
    //     data_file: impl AsRef<Path>,
    //     valid_samples: &Samples,
    //     min_snp_sep: u32,
    //     rec_rate: f64,
    // ) -> (Genome, Matrix<u8>, Sites) {
    //     let mut sites = Sites::new();
    //     let mut gbuilder = GenomeBuilder::new();

    //     let mut line = String::with_capacity(100000);
    //     let mut f = std::fs::File::open(data_file.as_ref())
    //         .map(BufReader::new)
    //         .unwrap();

    //     // get all sample names
    //     f.read_line(&mut line).unwrap();
    //     let samples: Vec<_> = line
    //         .trim()
    //         .split('\t')
    //         .skip(2)
    //         .map(|x| x.to_owned())
    //         .collect();
    //     line.clear();

    //     let n_valid_samples = samples
    //         .iter()
    //         .filter(|s| valid_samples.m().contains_key(*s))
    //         .count();

    //     let mut geno1 = MatrixBuilder::<u8>::new(n_valid_samples);
    //     let mut last_chrname = String::new();
    //     let mut last_pos = 0;
    //     while f.read_line(&mut line).unwrap() != 0 {
    //         let mut fields = line.trim().split("\t");
    //         let chrname = fields.next().unwrap();
    //         let pos: u32 = fields.next().unwrap().parse().unwrap();

    //         // println!("pos: {pos}, last_pos: {last_pos}");
    //         if (chrname == last_chrname) && (last_pos + min_snp_sep > pos) {
    //             line.clear();
    //             continue;
    //         } else {
    //             last_chrname.clear();
    //             last_chrname.push_str(chrname);
    //             last_pos = pos;
    //         }

    //         let (_, gw_pos) = gbuilder.encode_pos(chrname, pos);

    //         sites.add(gw_pos, gw_pos as f64 * rec_rate * 100.0);

    //         let mut cnt = 0;
    //         fields.enumerate().for_each(|(i, field)| {
    //             if valid_samples.m().contains_key(&samples[i]) {
    //                 let allel = match field.parse::<i8>().unwrap() {
    //                     -1 => None,
    //                     x => Some(x as u8),
    //                 };
    //                 geno1.push(allel);
    //                 cnt += 1;
    //             }
    //         });
    //         assert_eq!(cnt, n_valid_samples);
    //         line.clear();
    //     }

    //     let genome = gbuilder.finish();
    //     let geno1 = geno1.finish();
    //     sites.finish(&genome);

    //     (genome, geno1, sites)
    // }
    fn read_freq_file(
        freq_file: impl AsRef<Path>,
        args: &Arguments,
        genome: &Genome,
        sites: &Sites,
    ) -> Result<Matrix<f64>, Error> {
        let mut f = std::fs::File::open(freq_file.as_ref())
            .map(BufReader::new)
            .map_err(|source| Error::Io {
                source,
                file: Some(freq_file.as_ref().to_string_lossy().to_string()),
            })?;

        let mut freq = MatrixBuilder::<f64>::new(args.max_all as usize);
        let mut v = Vec::new();

        let mut line = String::with_capacity(100000);
        line.clear();

        let mut last_chrname = String::new();
        let mut last_pos = 0;
        while f
            .read_line(&mut line)
            .map_err(|_| ParseLineError::ReadLineError)?
            != 0
        {
            let mut fields = line.trim().split("\t");
            let chrname = fields
                .next()
                .ok_or(ParseLineError::ReadColumnError("chrname"))?;
            let pos: u32 = fields
                .next()
                .ok_or(ParseLineError::ReadColumnError("pos"))?
                .parse()
                .map_err(|_| ParseLineError::ParseColumnError("pos"))?;
            if (chrname == last_chrname) && (last_pos + args.min_snp_sep > pos) {
                line.clear();
                continue;
            } else {
                last_chrname.clear();
                last_chrname.push_str(chrname);
                last_pos = pos;
            }

            let gw_pos = genome.to_gw_pos(chrname, pos);
            assert!(sites.has_gw_pos(gw_pos));

            v.clear();
            v.resize(args.max_all as usize, 0.0);
            fields
                .enumerate()
                .try_for_each(|(i, field)| -> Result<(), Error> {
                    let af = field
                        .parse::<f64>()
                        .map_err(|_| ParseLineError::ParseColumnError("freq"))?;
                    v[i] = af;
                    Ok(())
                })?;
            for af in v.iter() {
                freq.push(Some(*af));
            }
            line.clear();
        }

        let freq = freq.finish();
        Ok(freq)
    }
    pub fn infer_freq_from_data(
        geno: &Matrix<u8>,
        sites: &Sites,
        args: &Arguments,
    ) -> Result<Matrix<f64>, Error> {
        // assert gentoeyps is still site oriented
        assert_eq!(geno.get_nrows(), sites.get_pos_slice().len());

        let mut freq = MatrixBuilder::<f64>::new(args.max_all as usize);
        let mut cnts = vec![0u32; args.max_all as usize];
        let mut total;

        for row in 0..geno.get_nrows() {
            cnts.clear();
            cnts.resize(args.max_all as usize, 0);
            total = 0;
            for allele in geno.get_row_iter(row).flatten() {
                total += 1;
                cnts[allele as usize] += 1;
            }
            if total == 0 {
                return Err(Error::FreqInferenceZeroNonMissGenotype);
            }

            for each in cnts.iter() {
                let af = Some(*each as f64 / total as f64);
                freq.push(af);
            }
        }

        Ok(freq.finish())
    }

    /// Infer allele frequencies from a subset of the samples of a genotype
    /// matrix.
    ///
    /// `cols` holds the column indices of the genotype matrix, that is, the
    /// samples the allele frequencies are calculated from. All samples of the
    /// genotype matrix are still analyzed by the HMM.
    pub fn infer_freq_from_data_subset(
        geno: &Matrix<u8>,
        sites: &Sites,
        genome: &Genome,
        args: &Arguments,
        cols: &[usize],
        freq_sample_file: &str,
    ) -> Result<Matrix<f64>, Error> {
        // assert gentoeyps is still site oriented
        assert_eq!(geno.get_nrows(), sites.get_pos_slice().len());

        let mut freq = MatrixBuilder::<f64>::new(args.max_all as usize);
        let mut cnts = vec![0u32; args.max_all as usize];
        let mut total;

        for row in 0..geno.get_nrows() {
            cnts.clear();
            cnts.resize(args.max_all as usize, 0);
            total = 0;
            let row_slice = geno.get_row_raw_slice(row);
            for col in cols.iter() {
                if let Some(allele) = row_slice[*col].as_option() {
                    total += 1;
                    cnts[allele as usize] += 1;
                }
            }
            if total == 0 {
                let (_chrid, chrname, pos) = genome.to_chr_pos(sites.get_pos_slice()[row]);
                return Err(Error::FreqInferenceZeroNonMissGenotypeInSubset {
                    file: freq_sample_file.to_owned(),
                    chrname: chrname.to_owned(),
                    pos,
                });
            }

            for each in cnts.iter() {
                let af = Some(*each as f64 / total as f64);
                freq.push(af);
            }
        }

        Ok(freq.finish())
    }

    /// Read a file of sample ids and map them to columns of a genotype matrix
    ///
    /// `col_offset` is the index, within all valid samples, of the first sample
    /// of the population the genotype matrix belongs to; `ncols` is the number
    /// of samples of that population. Ids that are not among the valid samples,
    /// for instance samples removed by the bcf filtering or by `-b/--bad-file`,
    /// are ignored with a warning.
    fn get_freq_sample_cols(
        freq_sample_file: &str,
        valid_samples: &Samples,
        col_offset: usize,
        ncols: usize,
    ) -> Result<Vec<usize>, Error> {
        let s = std::fs::read_to_string(freq_sample_file).map_err(|source| Error::Io {
            source,
            file: Some(freq_sample_file.to_owned()),
        })?;

        let mut cols = vec![];
        let mut not_found = vec![];
        for name in s.lines().map(|line| line.trim()).filter(|l| !l.is_empty()) {
            match valid_samples.m().get(name) {
                None => not_found.push(name),
                Some(idx) => {
                    let idx = *idx as usize;
                    if (idx < col_offset) || (idx >= col_offset + ncols) {
                        return Err(Error::FreqSampleWrongPopulation {
                            sample: name.to_owned(),
                            file: freq_sample_file.to_owned(),
                        });
                    }
                    cols.push(idx - col_offset);
                }
            }
        }
        cols.sort_unstable();
        cols.dedup();

        if !not_found.is_empty() {
            eprintln!(
                "WARN: {} of the {} sample ids listed in {} are not among the analyzed \
                samples and are ignored in the allele frequency calculation, including {:?}",
                not_found.len(),
                not_found.len() + cols.len(),
                freq_sample_file,
                &not_found[..not_found.len().min(5)],
            );
        }
        if cols.is_empty() {
            return Err(Error::FreqSamplesNotFound(freq_sample_file.to_owned()));
        }
        eprintln!(
            "allele frequencies are calculated from {} of the {} samples of this population",
            cols.len(),
            ncols
        );

        Ok(cols)
    }

    /// Restrict the genotypes of two populations read from two bcf/bin files to
    /// the sites they share
    ///
    /// Each file is filtered on its own, so the two sites sets are generally not
    /// identical, unlike for the text input format where identical sites are
    /// required. The shared sites are kept in the order of the first population
    /// and are thinned by the `--min-snp-sep` rule.
    fn intersect_sites_of_two_populations(
        geno1: Matrix<u8>,
        sitesinfo1: SiteInfoRaw,
        geno2: Matrix<u8>,
        sitesinfo2: SiteInfoRaw,
        min_snp_sep: u32,
    ) -> Result<(Matrix<u8>, Matrix<u8>, SiteInfoRaw), Error> {
        let nsites1 = sitesinfo1.get_chr_pos_vec().len();
        let nsites2 = sitesinfo2.get_chr_pos_vec().len();

        // map (chromosome name, position) of the second population to its row
        // in the second genotype matrix
        let mut sites2_map = HashMap::<(&str, u32), usize>::with_capacity(nsites2);
        for (irow, (chr_idx, pos)) in sitesinfo2
            .get_chr_idx_vec()
            .iter()
            .zip(sitesinfo2.get_chr_pos_vec().iter())
            .enumerate()
        {
            let chrname = sitesinfo2.get_chrname_vec()[*chr_idx].as_str();
            sites2_map.entry((chrname, *pos)).or_insert(irow);
        }

        let mut rows1 = vec![];
        let mut rows2 = vec![];
        let mut chr_idx_vec = vec![];
        let mut chr_pos_vec = vec![];
        let mut last_chrname = "";
        let mut last_pos = 0u32;
        for (irow, (chr_idx, pos)) in sitesinfo1
            .get_chr_idx_vec()
            .iter()
            .zip(sitesinfo1.get_chr_pos_vec().iter())
            .enumerate()
        {
            let chrname = sitesinfo1.get_chrname_vec()[*chr_idx].as_str();
            let irow2 = match sites2_map.get(&(chrname, *pos)) {
                Some(irow2) => *irow2,
                None => continue,
            };
            if (chrname == last_chrname) && (last_pos + min_snp_sep > *pos) {
                continue;
            }
            last_chrname = chrname;
            last_pos = *pos;
            rows1.push(irow);
            rows2.push(irow2);
            chr_idx_vec.push(*chr_idx);
            chr_pos_vec.push(*pos);
        }

        eprintln!(
            "two populations: n site = {} (population 1), {} (population 2), \
            {} shared sites are used",
            nsites1,
            nsites2,
            rows1.len()
        );
        if rows1.is_empty() {
            return Err(Error::NoSharedSitesBetweenPopulations);
        }

        let new_geno1 = Self::subset_geno_rows(&geno1, &rows1);
        let new_geno2 = Self::subset_geno_rows(&geno2, &rows2);

        let sitesinfo = SiteInfoRaw::from_parts(
            sitesinfo1.get_chrname_vec().to_vec(),
            sitesinfo1.get_chrname_map().clone(),
            chr_pos_vec,
            chr_idx_vec,
        );

        Ok((new_geno1, new_geno2, sitesinfo))
    }

    /// Keep only the given rows/sites of a site-oriented genotype matrix
    fn subset_geno_rows(geno: &Matrix<u8>, rows: &[usize]) -> Matrix<u8> {
        let ncols = geno.get_ncols();
        let mut data = Vec::<u8>::with_capacity(rows.len() * ncols);
        for irow in rows.iter() {
            data.extend_from_slice(geno.get_row_raw_slice(*irow));
        }
        Matrix::<u8>::from_shape_vec(rows.len(), ncols, data)
    }

    /// Keep only the given rows/sites of the site information
    fn subset_siteinfo_rows(sitesinfo: &SiteInfoRaw, rows: &[usize]) -> SiteInfoRaw {
        SiteInfoRaw::from_parts(
            sitesinfo.get_chrname_vec().to_vec(),
            sitesinfo.get_chrname_map().clone(),
            rows.iter()
                .map(|irow| sitesinfo.get_chr_pos_vec()[*irow])
                .collect(),
            rows.iter()
                .map(|irow| sitesinfo.get_chr_idx_vec()[*irow])
                .collect(),
        )
    }

    /// Remove the sites that have no genotype at all among the samples the
    /// allele frequencies are calculated from
    ///
    /// Restricting the frequency calculation to a subset of the samples, see
    /// `--freq-samples1`/`--freq-samples2`, can leave sites for which that
    /// subset carries no genotype and hence no allele frequency. Such sites
    /// cannot be scored by the HMM and are removed from the analysis.
    fn drop_sites_without_freq_samples(
        geno1: Matrix<u8>,
        geno2: Option<Matrix<u8>>,
        sitesinfo: SiteInfoRaw,
        freq_cols1: Option<&[usize]>,
        freq_cols2: Option<&[usize]>,
    ) -> (Matrix<u8>, Option<Matrix<u8>>, SiteInfoRaw) {
        if freq_cols1.is_none() && freq_cols2.is_none() {
            return (geno1, geno2, sitesinfo);
        }

        let has_genotype = |geno: &Matrix<u8>, cols: Option<&[usize]>, irow: usize| -> bool {
            match cols {
                None => true,
                Some(cols) => {
                    let row = geno.get_row_raw_slice(irow);
                    cols.iter().any(|icol| row[*icol].is_some())
                }
            }
        };

        let nsites = geno1.get_nrows();
        let rows: Vec<usize> = (0..nsites)
            .filter(|irow| {
                has_genotype(&geno1, freq_cols1, *irow)
                    && geno2
                        .as_ref()
                        .map_or(true, |geno2| has_genotype(geno2, freq_cols2, *irow))
            })
            .collect();

        if rows.len() == nsites {
            return (geno1, geno2, sitesinfo);
        }
        eprintln!(
            "WARN: {} of the {} sites have no genotype among the samples used for the \
            allele frequency calculation and are removed from the analysis",
            nsites - rows.len(),
            nsites
        );

        let new_geno1 = Self::subset_geno_rows(&geno1, &rows);
        let new_geno2 = geno2
            .as_ref()
            .map(|geno2| Self::subset_geno_rows(geno2, &rows));
        let new_sitesinfo = Self::subset_siteinfo_rows(&sitesinfo, &rows);

        (new_geno1, new_geno2, new_sitesinfo)
    }

    fn get_valid_pair_file(
        args: &Arguments,
        valid_samples: &Samples,
    ) -> Result<Vec<(u32, u32)>, Error> {
        let mut v = vec![];
        if let Some(good_file) = args.good_file.as_ref() {
            let mut buf = String::new();
            std::fs::File::open(good_file)
                .map(BufReader::new)
                .map_err(|e| Error::Io {
                    source: e,
                    file: Some(good_file.to_owned()),
                })?
                .read_to_string(&mut buf)
                .map_err(|e| Error::Io {
                    source: e,
                    file: Some(good_file.to_owned()),
                })?;
            let m = valid_samples.m();
            for line in buf.trim().split("\n") {
                let mut fields = line.split("\t");
                let s1 = fields.next().ok_or(ParseLineError::ReadColumnError("s1"))?;
                let s2 = fields.next().ok_or(ParseLineError::ReadColumnError("s2"))?;
                if m.contains_key(s1) && m.contains_key(s2) {
                    v.push((m[s1], m[s2]));
                }
            }
        } else {
            // two populatoin
            if valid_samples.pop2_nsam() > 0 {
                for i in 0..valid_samples.pop1_nsam() {
                    for j in valid_samples.pop1_nsam()..valid_samples.v().len() as u32 {
                        v.push((i, j))
                    }
                }
            }
            // one population
            else {
                let s1 = valid_samples.v();
                let n = s1.len() as u32;
                for i in 0..(n - 1) {
                    for j in (i + 1)..n {
                        v.push((i, j))
                    }
                }
            }
        }
        Ok(v)
    }

    /// get number of different alleles for all valide sites
    fn get_nall(freq1: &Matrix<f64>, freq2: Option<&Matrix<f64>>) -> Result<Vec<u8>, Error> {
        let mut v = vec![];
        // println!("{:?}", freq1.get_row_raw_slice(0));
        // println!("{:?}", freq2.unwrap().get_row_raw_slice(0));
        for i in 0..freq1.get_nrows() {
            let mut it1 =
                freq1
                    .get_row_iter(i)
                    .enumerate()
                    .map(|(iallele, x)| -> Result<bool, Error> {
                        Ok(x.ok_or(Error::MissingFrequency {
                            irow: i,
                            iallele,
                            value: x,
                        })? > 0.0)
                    });
            let n = match freq2 {
                Some(geno2) => {
                    let it2 = geno2.get_row_iter(i).enumerate().map(
                        |(iallele, x)| -> Result<bool, Error> {
                            Ok(x.ok_or(Error::MissingFrequency {
                                irow: i,
                                iallele,
                                value: x,
                            })? > 0.0)
                        },
                    );
                    it1.zip(it2)
                        .try_fold(0usize, |acc, (x, y)| -> Result<usize, Error> {
                            Ok(acc + (x? || y?) as usize)
                        })?
                }
                None => it1.try_fold(0usize, |acc, x| -> Result<usize, Error> {
                    Ok(acc + x? as usize)
                })?,
            };
            v.push(n as u8);
        }
        Ok(v)
    }

    /// Get the major allele for each site
    ///
    /// When one site has zero alleles, `Err(Error)` will be returned.
    fn get_major_all(
        freq1: &Matrix<f64>,
        freq2: Option<&Matrix<f64>>,
        nsam1_valid: usize,
        nsam2_valid: Option<usize>,
    ) -> Result<Vec<u8>, Error> {
        let mut v = vec![];
        let freq2 = freq2.unwrap_or(freq1);
        let nsam2_valid = nsam2_valid.unwrap_or(nsam1_valid) as f64;
        let nsam1_valid = nsam1_valid as f64;

        for ipos in 0..freq1.get_nrows() {
            let it = freq1.get_row_iter(ipos).map(|x| x.unwrap_or(0.0));
            let it2 = freq2.get_row_iter(ipos).map(|x| x.unwrap_or(0.0));

            let major_all = it
                .zip(it2)
                .map(|(a, b)| a * nsam1_valid + b * nsam2_valid)
                .enumerate()
                .try_fold(
                    (0, f64::MIN),
                    |(max_idx, max_val), (this_idx, this_val)| -> Result<(usize, f64), Error> {
                        if let std::cmp::Ordering::Greater = this_val
                            .partial_cmp(&max_val)
                            .ok_or(Error::PartialCmpIsNone("get_major_all"))?
                        {
                            Ok((this_idx, this_val))
                        } else {
                            Ok((max_idx, max_val))
                        }
                    },
                )?
                .0;

            // let major_all = it
            //     .zip(it2)
            //     .map(|(a, b)| a * nsam1_valid + b * nsam2_valid)
            //     .enumerate()
            //     .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            //     .ok_or(Error::EmptyIterator)?
            //     .0;
            v.push(major_all as u8);
            // if major_all > 1 {
            //     println!("{major_all}");
            // }
        }
        // println!("freq1.shape={},{}", freq1.get_nrows(), freq1.get_ncols());
        // println!("freq2.shape={},{}", freq2.get_nrows(), freq2.get_ncols());
        // println!("nsam1_valid={}, nsam2_valid={}", nsam1_valid, nsam2_valid);
        // println!("v.len={}", v.len());
        Ok(v)
    }
}

pub struct OutputFiles {
    pub seg_file: Arc<Mutex<BufWriter<File>>>,
    pub frac_file: Arc<Mutex<BufWriter<File>>>,
}

impl OutputFiles {
    pub fn new_from_args(
        args: &Arguments,
        buffer_size_segments: Option<usize>,
        buffer_size_frac: Option<usize>,
    ) -> Result<Self, Error> {
        let prefix = match args.output.as_ref() {
            Some(output) => output,
            None => &args.data_file1,
        };
        let seg_fn = format!("{prefix}.hmm.txt");
        let frac_fn = format!("{prefix}.hmm_fract.txt");
        use std::io::Write;

        let mut seg_file = match buffer_size_segments {
            Some(bfsz) => {
                std::fs::File::create(&seg_fn).map(|seg_fn| BufWriter::with_capacity(bfsz, seg_fn))
            }
            None => std::fs::File::create(&seg_fn).map(BufWriter::new),
        }
        .map_err(|e| Error::Io {
            source: e,
            file: Some(seg_fn.clone()),
        })?;

        let mut frac_file = match buffer_size_frac {
            Some(bfsz) => std::fs::File::create(&frac_fn)
                .map(|frac_fn| BufWriter::with_capacity(bfsz, frac_fn)),
            None => std::fs::File::create(&frac_fn).map(BufWriter::new),
        }
        .map_err(|e| Error::Io {
            source: e,
            file: Some(seg_fn.clone()),
        })?;

        // write header:
        writeln!(
            &mut frac_file,
            "sample1\tsample2\tN_informative_sites\tdiscordance\tlog_p\tN_fit_iteration\tN_generation\tN_state_transition\tseq_shared_best_traj\tfract_sites_IBD\tfract_vit_sites_IBD"
        )
        .map_err(|e| Error::Io {
            source: e,
            file: Some(frac_fn.clone()),
        })?;

        writeln!(
            &mut seg_file,
            "sample1\tsample2\tchr\tstart\tend\tdifferent\tNsnp"
        )
        .map_err(|e| Error::Io {
            source: e,
            file: Some(seg_fn.clone()),
        })?;

        Ok(Self {
            seg_file: Arc::new(Mutex::new(seg_file)),
            frac_file: Arc::new(Mutex::new(frac_file)),
        })
    }
}

#[test]
fn read_inputdata() {
    let mut args = Arguments::new_for_test();
    args.freq_file1 = None;
    let _input = InputData::from_args(&args);
}

#[test]
fn read_inputdata_freq_calc_error() {
    let mut args = Arguments::new_for_test();
    // set freq_file1 and data files
    args.freq_file1 = None;
    args.freq_file2 = None;
    args.data_file1 = "_tmp_hmm_input.txt".to_string();
    args.data_file2 = None;

    // write data file
    std::fs::write(
        &args.data_file1,
        concat!(
            "chr\tpos\ts1\ts2\ts3\ts4\n",
            "1\t100\t1\t0\t1\t0\n",
            "1\t200\t-1\t-1\t-1\t-1\n",
            "1\t300\t1\t0\t1\t0\n"
        ),
    )
    .unwrap();

    assert!(matches!(
        InputData::from_args(&args),
        Err(Error::FreqInferenceZeroNonMissGenotype)
    ));
}
#[test]
fn read_inputdata_freq_samples_subset() {
    let mut args = Arguments::new_for_test();
    args.freq_file1 = None;
    args.freq_file2 = None;
    args.data_file2 = None;

    // frequencies calculated from all samples
    let all = InputData::from_args(&args).unwrap();

    // frequencies calculated from the first 5 samples only
    let freq_sample_file = "_tmp_freq_samples1.txt";
    std::fs::write(freq_sample_file, all.samples.v()[..5].join("\n")).unwrap();
    args.freq_samples1 = Some(freq_sample_file.to_owned());
    let subset = InputData::from_args(&args).unwrap();

    // the same samples are analyzed (the genotype matrix is sample-oriented)
    assert_eq!(subset.geno.get_nrows(), all.geno.get_nrows());
    // and so are the same sites, except those without any genotype among the
    // samples the frequencies are calculated from
    assert!(subset.geno.get_ncols() <= all.geno.get_ncols());
    assert_eq!(subset.freq1.get_nrows(), subset.geno.get_ncols());
    assert_eq!(subset.freq1.get_nrows(), subset.sites.get_pos_slice().len());
    // but the allele frequencies are not the same
    assert!(subset.freq1.as_slice() != all.freq1.as_slice());

    // a file without any known sample id is an error
    std::fs::write(freq_sample_file, "not_a_sample_id").unwrap();
    assert!(matches!(
        InputData::from_args(&args),
        Err(Error::FreqSamplesNotFound(_))
    ));

    std::fs::remove_file(freq_sample_file).unwrap();
}

#[test]
fn read_inputdata_two_bcf_populations() {
    let mut args = Arguments::new_for_test_bcf();
    let bcf_filter_args =
        BcfFilterArgs::new_from_toml_file(args.bcf_filter_config.as_ref().unwrap()).unwrap();

    // split the samples of the test bcf into two populations and write the
    // genotypes of each population to its own binary file
    let all =
        BcfGenotype::new_from_processing_bcf(&args.bcf_read_mode, &bcf_filter_args, &args.data_file1)
            .unwrap();
    let samples = all.get_samples().to_vec();
    assert!(samples.len() >= 4);
    let (pop1, pop2) = samples.split_at(samples.len() / 2);

    let mut bin_files = vec![];
    for (ipop, pop) in [pop1, pop2].iter().enumerate() {
        let sample_file = format!("_tmp_two_pop{ipop}_samples.txt");
        let bin_file = format!("_tmp_two_pop{ipop}.bin");
        std::fs::write(&sample_file, pop.join("\n")).unwrap();

        let mut bcf_filter_args = bcf_filter_args.clone();
        bcf_filter_args.target_samples = Some(sample_file.clone().into());
        BcfGenotype::new_from_processing_bcf(
            &args.bcf_read_mode,
            &bcf_filter_args,
            &args.data_file1,
        )
        .unwrap()
        .save_to_file(&bin_file)
        .unwrap();

        std::fs::remove_file(&sample_file).unwrap();
        bin_files.push(bin_file);
    }

    args.from_bcf = false;
    args.from_bin = true;
    args.bcf_filter_config = None;
    args.data_file1 = bin_files[0].clone();
    args.data_file2 = Some(bin_files[1].clone());

    let input = InputData::from_args(&args).unwrap();

    // both populations are read and the second one gets its own frequencies
    assert!(input.samples.pop1_nsam() > 0);
    assert!(input.samples.pop2_nsam() > 0);
    assert!(input.freq2.is_some());
    // the genotype matrix is sample-oriented and holds both populations
    assert_eq!(
        input.geno.get_nrows() as u32,
        input.samples.pop1_nsam() + input.samples.pop2_nsam()
    );
    // sites are shared by the two populations
    let nsites = input.sites.get_pos_slice().len();
    assert!(nsites > 0);
    assert_eq!(input.geno.get_ncols(), nsites);
    assert_eq!(input.freq1.get_nrows(), nsites);
    assert_eq!(input.freq2.as_ref().unwrap().get_nrows(), nsites);
    // cross population pairs only
    assert!(!input.get_chunk_pairs().is_empty());

    for bin_file in bin_files.iter() {
        std::fs::remove_file(bin_file).unwrap();
    }
}

pub struct FracRecord<'a> {
    pub sample1: &'a str,
    pub sample2: &'a str,
    pub sum: usize,
    pub discord: f64,
    pub max_phi: f64,
    pub iter: usize,
    pub k_rec: f64,
    pub ntrans: usize,
    pub seq_ibd_ratio: f64,
    pub count_ibd_fb_ratio: f64,
    pub count_ibd_vit_ratio: f64,
}

pub struct SegRecord<'a> {
    pub sample1: &'a str,
    pub sample2: &'a str,
    pub chrname: &'a str,
    pub start_pos: u32,
    pub end_pos: u32,
    pub ibd: u8,
    pub n_snp: usize,
}

pub struct OutputBuffer<'a> {
    seg_file: Arc<Mutex<BufWriter<File>>>,
    frac_file: Arc<Mutex<BufWriter<File>>>,
    segs: SmallVec<[SegRecord<'a>; 5]>,
    fracs: SmallVec<[FracRecord<'a>; 1]>,
}

impl<'a> OutputBuffer<'a> {
    pub fn new(out: &OutputFiles, segs_capacity: usize, fracs_capacity: usize) -> Self {
        Self {
            seg_file: Arc::clone(&out.seg_file),
            frac_file: Arc::clone(&out.frac_file),
            segs: SmallVec::<[SegRecord<'a>; 5]>::with_capacity(segs_capacity),
            fracs: SmallVec::<[FracRecord<'a>; 1]>::with_capacity(fracs_capacity),
        }
    }

    pub fn add_seg(&mut self, seg: SegRecord<'a>) -> Result<(), Error> {
        if self.segs.len() == self.segs.capacity() {
            self.flush_segs()?;
        }
        self.segs.push(seg);
        Ok(())
    }

    pub fn add_frac(&mut self, frac: FracRecord<'a>) -> Result<(), Error> {
        if self.fracs.len() == self.fracs.len() {
            self.flush_frac()?;
        }
        self.fracs.push(frac);
        Ok(())
    }

    pub fn flush_segs(&mut self) -> Result<(), Error> {
        use std::io::Write;
        let mut file = self
            .seg_file
            .lock()
            .map_err(|_| Error::LockError("flush_segs"))?;
        for seg in self.segs.iter() {
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                seg.sample1,
                seg.sample2,
                seg.chrname,
                seg.start_pos,
                seg.end_pos,
                seg.ibd,
                seg.n_snp,
            )
            .unwrap();
        }
        self.segs.clear();
        Ok(())
    }

    pub fn flush_frac(&mut self) -> Result<(), Error> {
        use std::io::Write;
        let mut file = self
            .frac_file
            .lock()
            .map_err(|_| Error::LockError("flush_frac"))?;
        for frac in self.fracs.iter() {
            writeln!(
                file,
                "{}\t{}\t{}\t{:.4}\t{:0.5e}\t{}\t{:.3}\t{}\t{:.5}\t{:.5}\t{:.5}",
                frac.sample1,
                frac.sample2,
                frac.sum,
                frac.discord,
                frac.max_phi,
                frac.iter,
                frac.k_rec,
                frac.ntrans,
                frac.seq_ibd_ratio,
                frac.count_ibd_fb_ratio,
                frac.count_ibd_vit_ratio,
            )
            .unwrap();
        }
        self.fracs.clear();
        Ok(())
    }
}
