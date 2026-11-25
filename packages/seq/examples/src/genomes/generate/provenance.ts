/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { GenomesFiles } from './index.js';

// The Genomes dataset's generated files as `sources.json` describes them (licence `LicenseRef-seqx-generated`): what each
// is and how it is derived. The record script's `genomes` case writes a row a file from these; the two files folded in
// from the formats' generators (`generated-genes.bed`, `generated-genotypes.vcf`) keep their own treatment (compared
// with their generators, no row).

type Described = Exclude<keyof typeof GenomesFiles, 'genesBed' | 'genotypes'>;

const From = '(src/genomes/generate/, seeded on core’s PCG and hash3, a stream per seed, concern and index)';

export const GenomesProvenance: Readonly<Record<Described, { readonly entry: string; readonly derived: string }>> = {
  sizesA: {
    entry: 'toy A’s sequences and lengths as UCSC chrom.sizes: chr1 to chr3, chrM (cyclic) and HLA-A*01:01:01:01',
    derived: `the toy assembly’s lengths ${From}`,
  },
  sizesB: {
    entry: 'toy B’s sequences (unprefixed: 1, 2, 3, MT) and lengths as chrom.sizes',
    derived: `toy B made from toy A by the alignment truth (an inversion, a translocation, a gap, small indels) ${From}`,
  },
  sizesA2: {
    entry: 'toy A′ (toy A revised: a duplication, an insertion, a deletion) as chrom.sizes',
    derived: `toy A’s lengths with the revision’s events ${From}`,
  },
  fasta: {
    entry: 'toy A’s chrM (16,569) and HLA-A*01:01:01:01 (3,503) as FASTA, 60 letters a line',
    derived: `letters a 4 kb tile at a time from a stream per sequence and tile ${From}`,
  },
  aliases: {
    entry:
      'toy A’s alias table in the seq format (coverage only): each chromosome’s bare number and an NC_-style accession, chrM its accession (no MT)',
    derived: `the toy assembly’s names ${From}`,
  },
  alignment: {
    entry:
      'toy A against toy B: gapless blocks as sided pairs in the seq format (side a toy A, side b toy B; row the alignment)',
    derived: `the alignment truth: blocks of 3 to 9 kb, a small indel between two, an 80 kb inversion, a 60 kb translocation, chr1’s run of N unaligned ${From}`,
  },
  dots: {
    entry:
      'the dot plot of toy A against toy B as sparse grid tiles in the seq format (MessagePack), cells of 1,024 and 4,096 positions',
    derived: `the alignment truth’s positions counted in cells ${From}`,
  },
  liftover: {
    entry: 'toy A to toy A′: gapless runs as sided pairs in the seq format (side a A, side b A′, named alike)',
    derived: `the revision’s runs: 10 kb of chr1 held twice, 5 kb inserted in chr2, 4 kb of chr3 lost; at level 1 and as bin pairs of 16,384 ${From}`,
  },
  gff3: {
    entry: '150 genes on toy A in GFF3 (gene, mRNA, exon, CDS), bgzipped',
    derived: `genes of 3 to 12 exons and 1 to 4 transcripts, both strands, named TOY0001 to TOY0150 ${From}`,
  },
  gff3Index: {
    entry: 'the GFF3’s tabix index (TBI)',
    derived: `written as tabix writes it over the bgzipped GFF3 ${From}`,
  },
  gtf: {
    entry: 'the same 150 genes in GTF (gene, transcript, exon, CDS), bgzipped',
    derived: `the GFF3’s model written as GTF ${From}`,
  },
  gtfIndex: {
    entry: 'the GTF’s tabix index (TBI)',
    derived: `written as tabix writes it over the bgzipped GTF ${From}`,
  },
  signal: {
    entry: 'a TPM-like signal over toy A as a bedGraph of 100 bp bins (0 to 10^4), bgzipped',
    derived: `waves in log space and noise, no line over chr1’s run of N and 10 kb of chr2 (missing values) ${From}`,
  },
  signalIndex: {
    entry: 'the bedGraph’s tabix index (TBI)',
    derived: `written as tabix writes it over the bgzipped bedGraph ${From}`,
  },
  vcf: {
    entry: '40 samples’ genotypes at 2,000 sites of toy A as VCF 4.2 (GT), bgzipped',
    derived: `SNVs at toy A’s letters, chrM’s calls haploid, phased and unphased lines, partial and no-calls, a 600-base deletion over 12 SNVs carried by six samples ${From}`,
  },
  vcfIndex: {
    entry: 'the VCF’s tabix index (TBI)',
    derived: `written as tabix writes it over the bgzipped VCF ${From}`,
  },
  methylation: {
    entry:
      '20 samples’ methylation at 250 CpGs of toy A’s chr1 as a matrix in the seq format (MessagePack): items and 1,024-position bins',
    derived: `toy A’s own CG positions in a 50 kb stretch, an island near zero, 8 % missing (NaN) ${From}`,
  },
  geneMethylation: {
    entry:
      'the 150 genes’ methylation in 20 samples as a matrix in the seq format (MessagePack): items and 4,096-position bins',
    derived: `a level a gene, a quarter of the genes lower in the blood samples ${From}`,
  },
  contacts: {
    entry:
      'a symmetric Hi-C over toy A’s chr1 to chr3 as grid tiles in the seq format (MessagePack), cells of 4,096 and 16,384 positions',
    derived: `counts falling with the distance, domains along the diagonal, a hotspot between chr1 and chr3 ${From}`,
  },
};
