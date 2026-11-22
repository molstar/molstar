/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FilesViews, GenomeFormatsViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';

const view = (id: string) => GenomeFormatsViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Features/Genome formats/GenBank' };
export default meta;

/** The GenBank view and, under it, the GFF3 view of the same plasmid: two stories' roots in one page. */
function besideGff3(): HTMLElement {
  const root = document.createElement('div');
  root.style.cssText = 'display: grid; gap: var(--msui-space-6, 24px); background: var(--msui-bg);';
  const gff3 = FilesViews.find((v) => v.id === 'plasmid-gff3')!;
  root.append(exampleStory(view('genbank-plasmid'), { files: true }), exampleStory(gff3, { files: true }));
  return root;
}

/**
 * The plasmid’s generated GenBank record (pUC19c turned to an origin inside bla; `fixtures/plasmid/pUC19c-rotated.gb`)
 * read by the opt-in `genbank` runtime: its letters and its features from the one file, a track a feature key. The
 * record says it is circular; bla’s `complement(join(2312..2686,1..486))` and the primer site’s join over the origin
 * are one item each across it. Under it the GFF3 view of the same plasmid (Features/Files, letters and features from
 * two files): the two draw the same (a headless test compares their draw lists). The file’s records and its head are
 * printed under the view.
 */
export const Plasmid: StoryObj = { render: besideGff3 };

/**
 * A misread GenBank file: what parses is drawn; a remote reference, a site between two bases, a location past the
 * record’s length and a malformed one are each left out (the feature, never the file), an `a.b` location read as its
 * hull, more letters than LOCUS says (the letters win), a CONTIG assembly not read and a last record cut short; each
 * said once with its line numbers (the status line under the view).
 */
export const MisreadFile: StoryObj = {
  render: () => exampleStory(view('genbank-misread'), { files: true, sourceIssues: true }),
};
