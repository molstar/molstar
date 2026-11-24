/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GenomeFormatsViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { everyLane, findAndFocus } from './names.js';

const view = GenomeFormatsViews.find((v) => v.id === 'gtf-tp53')!;

const meta: Meta = { title: 'Features/Genome formats/GTF' };
export default meta;

/**
 * TP53’s region of GENCODE v50’s GTF (the recorded slice, gzip, read whole by the opt-in `gff` runtime): the same gene
 * models from GTF’s `transcript_id` and `gene_id` (its `transcript` lines give each extent), the thick part reaching
 * the stop codon GTF writes after the CDS, and one model a gene (`features: genes`). chr17’s length from an inline
 * chrom.sizes (a GTF says none).
 */
export const Tp53: StoryObj = { render: () => exampleStory(view, { files: true }) };

/**
 * A name typed: “tp53-201” found by the file’s name index (GTF’s `gene_name`, `gene_id`, `transcript_name`,
 * `transcript_id`, case-insensitive), printed and focused.
 */
export const FindByName: StoryObj = {
  render: () =>
    exampleStory(everyLane(view), {
      summary:
        'A name typed in the find box above the view is looked up in the GTF’s name index (`gene_name`, `gene_id`, `transcript_name`, `transcript_id`, case-insensitive): the first hit’s extent is focused, and what was searched, what was found and where the view went are printed under the box. Type another name (TP53, WRAP53-201, ENST00000269305) and press Find.',
      find: { title: 'Find a name', initial: 'tp53-201', run: (v, text) => findAndFocus(v, 'gtf', text) },
    }),
};
