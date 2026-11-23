/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ExampleView } from '../views.js';
import { namesAliased } from './by-region.js';
import * as compare from './compare.js';
import * as tracks from './tracks.js';
import { registry } from './toys.js';

// The Genomes dataset's views (examples.md, Genomes, The views), in the design's cut order reversed: a cut drops the
// last (the liftover, then the contact map). The genome ring and both genomes on one ring are the ring extension's.

const view = (id: string, title: string, summary: string, build: () => ReturnType<ExampleView['build']>) =>
  ({ dataset: 'genomes', id, title, summary, build, registry }) satisfies ExampleView;

export const GenomesViews: readonly ExampleView[] = [
  view(
    'synteny',
    'Toy A over toy B',
    'Two generated toy genomes, A (chr1, chr2, chr3, chrM, an HLA contig) above and B (named 1, 2, 3, MT) below, their alignment’s blocks as ribbons (sided pairs in the `seq` format, read by region from either end): A’s chr1 inversion crossed, A’s chr3 tail carried to B’s 2 (a translocation), chr1’s run of N no ribbon (a gap). Hover a ribbon: both ends light.',
    () => compare.synteny(),
  ),
  view(
    'dot-plot',
    'Toy A against toy B, a dot plot',
    'The same alignment as a dot plot: A across, B down, a grid following both sections, the aligned positions counted in cells of 1,024 and 4,096 positions (two levels of tiles) and the blocks on top as diagonals: the inversion an anti-diagonal, the translocation off the diagonal, the gap a break. Hover a cell: both positions; drag a rectangle: both ranges.',
    () => compare.dotPlot(),
  ),
  view(
    'signal',
    'Signal by region',
    'A TPM-like bedGraph over toy A (100 bp bins, 0 to 10^4) read by region through its tabix index at the coarsest level that fits a device column, three ways: its domain open (the data in view), open at its upper end only (from 0), and on a log scale. chr1’s run of N and 10 kb of chr2 have no value: gaps, never zero.',
    () => tracks.signal(),
  ),
  view(
    'samples',
    'Samples’ variants as rows',
    '40 samples’ genotypes on toy A’s chr1 around a 600-base deletion, read by region through the VCF’s tabix index: the sites (a row an ALT allele) and the matrix below, a row a sample coloured by genotype class (uncalled in the missing colour); the deletion one item solid across the 12 SNVs it covers in its six carriers’ rows. Hover a cell: the sample and the allele.',
    () => tracks.samples(),
  ),
  view(
    'methylation',
    'Methylation by region',
    'A methylation fraction per CpG and sample (20 samples, 250 CpGs of chr1, 8 % missing: left clear, never 0) with the island where it falls; below, the same matrix layer over the 150 genes, each gene’s mean a sample: a gene-by-sample heatmap over the genome.',
    () => tracks.methylation(),
  ),
  view(
    'transcripts',
    'Transcripts in three levels',
    'Toy A’s genes from a GFF3 read by region through its index: one model a gene, then its transcripts (gene → transcript → exon: exons as boxes, the CDS thick, a hairline with chevrons on the minus strand), and the GTF of the same model beside them. Hover a transcript: every exon lights.',
    () => tracks.transcripts(),
  ),
  view(
    'by-name',
    'Open on a gene by name',
    `A GFF3 read whole, its genes found by name (exact names first, then by prefix, over \`Name\` and \`ID\`): the story finds ${tracks.OpenedGene} and focuses its extent, known from the hit without reading a tile.`,
    () => tracks.byName(),
  ),
  view(
    'names',
    'Toy A by region, its names aliased',
    'Toy A’s letters by region from an in-memory adapter (no file: sequences, lengths and letters from the seed) with its alias table: marks written as `1`, `NC_900002.1` and `chr3` drawn on chr1 to chr3; one written `MT` refused (never guessed onto chrM). The HLA contig, `HLA-A*01:01:01:01`, read from a small FASTA.',
    () => namesAliased(),
  ),
  view(
    'contact-map',
    'Hi-C over toy A',
    'A symmetric contact map over toy A’s chr1 to chr3 read in tiles of two levels (the upper triangle served, the lower its mirror), on a log palette: domains along the diagonal, a hotspot between chr1 and chr3.',
    () => compare.contactMap(),
  ),
  view(
    'liftover',
    'Toy A lifted to A′',
    'Toy A’s genes lifted to its revision A′ through the map between them (gapless runs as sided pairs in the `seq` format, read by region): chr1’s duplicated 10 kb fanning to two places (its gene drawn twice), chr2 shifted by its 5 kb insertion, chr3’s lost 4 kb no ribbon. Hover a gene: its images light below.',
    () => compare.liftover(),
  ),
];
