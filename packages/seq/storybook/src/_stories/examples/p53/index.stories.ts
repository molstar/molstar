/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';

/** The views' stories: the title's last step and the export's name (Storybook's id: the sanitized title, `--`, the name). */
const Stories: Readonly<Record<string, readonly [title: string, name: string]>> = {
  tracks: ['Tracks', 'tracks'],
  'protein-on-gene': ['Protein on its gene', 'protein-on-its-gene'],
  'one-axis': ['Protein and gene on one axis', 'protein-and-gene-on-one-axis'],
  'genomic-follower': ['Genomic follower', 'genomic-follower'],
  'reversed-gene': ['Reversed gene', 'reversed-gene'],
  'gene-by-region': ['Gene by region', 'gene-by-region'],
  'edit-and-undo': ['Edit and undo', 'edit-and-undo'],
  'bookmarks-and-links': ['Bookmarks and links', 'bookmarks-and-links'],
  selection: ['Selection', 'selection'],
};

/** The examples built on the views: what to try in each. */
const Tries: readonly { readonly id: string; readonly title: string; readonly summary: string }[] = [
  {
    id: 'edit-and-undo',
    title: 'Edit and undo',
    summary: 'The tracks edited: add and hide a track, set a theme, undo and redo, drag a height to preview it.',
  },
  {
    id: 'bookmarks-and-links',
    title: 'Bookmarks and links',
    summary: 'The gene’s region with bookmarks and a text whose links move the view, on transcripts that arrive later.',
  },
  {
    id: 'selection',
    title: 'Selection',
    summary: 'Drag on the gene to select nucleotides; the protein shows them projected, in canvas and SVG.',
  },
];

const sanitize = (s: string) =>
  s
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-|-$/g, '');

function storyId(view: string): string {
  const [title, name] = Stories[view];
  return `${sanitize(`Seq* (next)/Examples/p53/${title}`)}--${name}`;
}

const meta: Meta = { title: 'Examples/p53/Index' };
export default meta;

/** The views of the p53 dataset, each a story of its own. */
export const Views: StoryObj = {
  render: () => {
    const root = document.createElement('div');
    root.style.cssText = 'font: 14px/1.5 system-ui, sans-serif; padding: 24px; max-width: 72ch;';
    const h = document.createElement('h2');
    h.textContent = 'p53';
    const intro = document.createElement('p');
    intro.textContent =
      'Human TP53: its protein with domains and conservation, and a toy gene it is placed on. Each view is one spec over this data.';
    const list = document.createElement('ul');
    for (const v of P53Views) {
      const li = document.createElement('li');
      const a = document.createElement('a');
      a.href = `?path=/story/${storyId(v.id)}`;
      a.target = '_top';
      a.textContent = v.title;
      li.append(a, ` ${v.summary}`);
      list.append(li);
    }
    const h2 = document.createElement('h3');
    h2.textContent = 'Editing, addressing and selection';
    const list2 = document.createElement('ul');
    for (const v of Tries) {
      const li = document.createElement('li');
      const a = document.createElement('a');
      a.href = `?path=/story/${storyId(v.id)}`;
      a.target = '_top';
      a.textContent = v.title;
      li.append(a, ` ${v.summary}`);
      list2.append(li);
    }
    root.append(h, intro, list, h2, list2);
    return root;
  },
};
