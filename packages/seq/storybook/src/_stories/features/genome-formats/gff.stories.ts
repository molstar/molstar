/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GenomeFormatsViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { findAndFocus } from './names.js';

const view = (id: string) => GenomeFormatsViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Features/Genome formats/GFF3' };
export default meta;

/**
 * TP53’s region of Ensembl GRCh38 release 116’s GFF3 (the recorded slice, gzip, read whole by the opt-in `gff`
 * runtime): its transcripts as gene models (the ranges layer’s `line: 'gene'`: exons as boxes, the CDS thick, the UTRs
 * thin, a hairline with chevrons pointing left for the minus strand), and from the same file with `features: genes` one
 * model a gene. The file’s head and the source nodes (the `features` param of each) are printed under the view.
 */
export const Tp53: StoryObj = { render: () => exampleStory(view('gff-tp53'), { files: true }) };

/**
 * A name typed: “TP53-2” found by the file’s name index (`findNames`, exact names first, then by prefix; GFF3’s `Name`,
 * `ID`, `Alias`), the hits printed under the view and the first focused (its transcript’s extent).
 */
export const FindByName: StoryObj = {
  render: () =>
    exampleStory(view('gff-tp53'), {
      aside: { title: 'Found by name', text: (v) => findAndFocus(v, 'gff3', 'TP53-2') },
    }),
};

/**
 * A misread GFF3: what parses is drawn as gene models; an end before its start, a position that is no number and a
 * short line are left out, an exon whose parent no line has and two transcripts each other’s parent are kept as items
 * of their own, a gene with a CDS and no transcript is read as one; each said once with its line number (the status
 * line under the view).
 */
export const MisreadFile: StoryObj = {
  render: () => exampleStory(view('gff-misread'), { files: true, sourceIssues: true }),
};
