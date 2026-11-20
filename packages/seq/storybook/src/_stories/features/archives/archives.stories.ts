/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ExampleViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { HostileArchives } from './archive-model.js';
import { archiveStory, refusalStory } from './archive-root.js';

const view = (dataset: string, id?: string) => ExampleViews.find((v) => v.dataset === dataset && (!id || v.id === id))!;

const meta: Meta = { title: 'Features/Archives' };
export default meta;

/**
 * The plasmid saved with its data (the construct’s `seq` file) as one zip, then opened again: the view under the
 * listing is drawn from the archive alone, its source read as an entry of it through the assets port, with the page
 * never asked for the file again.
 */
export const PlasmidSaved: StoryObj = {
  render: () => archiveStory(view('plasmid'), 'The plasmid view, saved with its data and opened from the archive.'),
};

/** p53 on its gene: the gene’s file is the entry the spec’s source names, under the path it had. */
export const P53OnItsGeneSaved: StoryObj = {
  render: () =>
    archiveStory(
      view('p53', 'protein-on-gene'),
      'p53 on its gene, saved with the gene’s file and opened from the archive.',
    ),
};

/**
 * Archives no reader opens (EC-91): an entry above the root and an absolute path are refused by name, a megabyte of
 * zeros in a few hundred bytes by its ratio, more entries than the bound by the count: each before anything is
 * inflated, with the issue the reader gives.
 */
export const Refused: StoryObj = { render: () => refusalStory(HostileArchives) };
