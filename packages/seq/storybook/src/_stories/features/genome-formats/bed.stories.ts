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

const meta: Meta = { title: 'Features/Genome formats/BED' };
export default meta;

/**
 * UCSC’s hg38 cytobands (the recorded file, gzip, read whole by the opt-in `bed` runtime): chr17’s bands drawn as bands
 * and as named ranges, the ideogram above showing the window. The band’s stain is the file’s fifth column, named by its
 * number (`fields: { stain: '4' }`); chr17’s length comes from an inline chrom.sizes. The file’s head and the source
 * nodes are printed under the view.
 */
export const Cytobands: StoryObj = story('bed-cytobands');

/** Generated BED12 genes: each exon a block, the blocks one item (hover one: the whole gene lights), strands drawn. */
export const Bed12Genes: StoryObj = story('bed-genes');

/**
 * A misread BED: the lines that parse are drawn; an end before its start and a position that is no number are left out,
 * a strand and a score that are not ones read as missing, spaces and an extra column said, each once with its line
 * number (the status line under the view).
 */
export const MisreadFile: StoryObj = story('bed-misread', true);
