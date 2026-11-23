/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GenomeFormatsViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';

const view = (id: string) => GenomeFormatsViews.find((v) => v.id === id)!;
const story = (id: string, sourceIssues = false): StoryObj => ({
  render: () => exampleStory(view(id), { files: true, sourceIssues }),
});

const meta: Meta = { title: 'Features/Genome formats/VCF' };
export default meta;

/**
 * TP53’s variants in NCBI ClinVar (GRCh38, the recorded slice: BGZF read by region at level 1 through the tabix index
 * the source names) by the opt-in `vcf` runtime: a row an ALT allele, a mark at its REF extent coloured by its
 * category (a run of marks merged under a pixel keeps its category’s colour, a mixed run is neutral). Hover an allele:
 * the card names it by ClinVar’s id, with its track and extent; click selects its REF extent, a double-click focuses it. The file’s header summary and the source nodes are printed under the view.
 */
export const Tp53ClinVar: StoryObj = story('vcf-clinvar');

/**
 * A generated file read whole: SNVs, multi-allelic sites, indels and MNVs, and two long deletions (a symbolic `<DEL>`
 * to its `END`, a 600-base sequence-resolved one), each one item drawn as a bar across the SNVs it covers.
 */
export const LongDeletion: StoryObj = story('vcf-long-deletion');

/**
 * Genotypes of 40 generated samples (read whole) by the `vcf` runtime’s matrix: a row a sample in the header’s order,
 * named in the label column (hover: `Sample S11 · 0/1 (heterozygous)`), a
 * column an ALT allele at its extent coloured by its genotype class (no copy, some, every copy, partly called; uncalled
 * in the missing colour); the last eight samples haploid; one long deletion one item across the SNVs it covers, solid
 * in its carriers’ rows (EC-220, EC-222). The sites are drawn above from the same rows; genotypes are read for the
 * matrix alone.
 */
export const Genotypes: StoryObj = story('vcf-genotypes');

/**
 * ClinVar’s chr17 at the fit: more level-1 tiles than a view reads, so none is asked: the track stays the band drawn
 * while data is waiting (no text, no issue in the status line) and nothing is read (no estimate drawn). Only TP53’s pages are recorded: zoom into 7.67–7.69 Mb to read them.
 */
export const ZoomIn: StoryObj = story('vcf-zoom-in');

/**
 * A misread VCF: what parses is drawn; a position that is no number, no REF and a short line are left out, a QUAL and
 * a depth that are not numbers read as missing, an undeclared INFO key, an `*` allele, a record past its contig
 * (clipped) and one beyond it (left out), and a header line without its Type said, each once with its line numbers
 * (the status line under the view).
 */
export const MisreadFile: StoryObj = story('vcf-misread', true);
