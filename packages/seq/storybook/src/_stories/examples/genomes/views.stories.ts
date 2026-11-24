/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GenomesViews } from '@molstar/seqx-examples';
import { OpenedGene } from '@molstar/seqx-examples/genomes/tracks';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';
import { findGeneAndFocus } from './find.js';

const view = (id: string) => GenomesViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/Genomes/Views' };
export default meta;

/**
 * Toy A over toy B (the dataset's generated genomes; B names its sequences unprefixed), their alignment's blocks as
 * ribbons, read by region from either end of a `seq` file of sided pairs: A's chr1 inversion crossed, A's chr3 tail on
 * B's 2 (a translocation), chr1's run of N no ribbon (a gap). The ribbons coloured by orientation, keyed in the
 * legend ("A to B": same, inverted; a ribbon's width its aligned length). Hover a ribbon: both ends light.
 */
export const Synteny: StoryObj = { render: () => exampleStory(view('synteny'), { files: true }) };

/**
 * The same alignment as a dot plot: A across, B down, two levels of grid tiles (the positions aligned in a cell) under
 * the blocks drawn as diagonals. Hover a cell: the crosshair and its two positions. Drag a rectangle: both ranges.
 */
export const DotPlot: StoryObj = { render: () => exampleStory(view('dot-plot')) };

/**
 * A bedGraph read by region through its tabix index, three ways: an open domain, `[0, null]`, a log scale (TPM 0.1 to
 * 10^4). Zoom in: a finer level is read, the coarser one drawn while it waits.
 */
export const Signal: StoryObj = { render: () => exampleStory(view('signal'), { files: true }) };

/**
 * 40 samples by a VCF read by region (its index named by the source): the sites, and the genotype matrix, the 600-base
 * deletion solid across the SNVs it covers in its carriers' rows. Hover a cell; click selects the allele.
 */
export const Samples: StoryObj = { render: () => exampleStory(view('samples'), { files: true }) };

/** CpG methylation (missing values clear, never 0) and the gene-by-sample heatmap of the same layer. */
export const Methylation: StoryObj = { render: () => exampleStory(view('methylation')) };

/** Genes, transcripts and the GTF of the same model, read by region: gene → transcript → exon. */
export const Transcripts: StoryObj = { render: () => exampleStory(view('transcripts'), { files: true }) };

/** A gene found by name in the GFF3 read whole: the first hit's extent focused. */
export const ByName: StoryObj = {
  render: () =>
    exampleStory(view('by-name'), {
      find: { title: 'Find a gene', initial: OpenedGene, run: (v, text) => findGeneAndFocus(v, text) },
    }),
};

/**
 * Toy A from the letters adapter (no file) with its alias table: marks named `1`, `NC_900002.1` and `chr3` drawn, one
 * named `MT` refused; the HLA contig from the small FASTA. Zoom in for the letters.
 */
export const NamesAliased: StoryObj = { render: () => exampleStory(view('names'), { sourceIssues: true }) };

/** A symmetric Hi-C over toy A's chr1 to chr3, two levels of tiles: domains, a hotspot between chr1 and chr3. */
export const ContactMap: StoryObj = { render: () => exampleStory(view('contact-map')) };

/**
 * A Manhattan-like plot: a toy association study's SNPs as dots at their height (-log10 p) over one 20 kb chromosome,
 * a dashed reference line at the genome-wide threshold (7.3) and a dotted one at the suggestive one (5), each labelled
 * with its value ("genome-wide · 7.3") at the right end. The strongest signal reaches about 6 and the domain is open, so the lines are part of it: the plot
 * shows that none passed the threshold.
 */
export const Manhattan: StoryObj = { render: () => exampleStory(view('manhattan')) };

/** Toy A's genes lifted to A′ through the map between them (a duplication drawn twice, a deletion, an insertion). */
export const Liftover: StoryObj = { render: () => exampleStory(view('liftover')) };
