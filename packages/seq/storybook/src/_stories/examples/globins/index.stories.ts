/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews, GlobinsViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { storyId } from '../../story-id.js';

// The Globins dataset's index (examples-and-stories): its views (each a story), the features shown on its data, and
// the sources with their licences and attribution (UniProt's is CC BY 4.0: shown here, and carried by the specs'
// source nodes).

const meta: Meta = { title: 'Examples/Globins/Index' };
export default meta;

/** Where each Globins view's story is: its title under `Seq* (next)` and its export. */
const Where: Readonly<Record<string, readonly [string, string]>> = {
  family: ['Examples/Globins/Rows', 'Family'],
  folded: ['Examples/Globins/Rows', 'Folded'],
  figure: ['Examples/Globins/Rows', 'Figure'],
  hemoglobin: ['Examples/Globins/Hemoglobin', 'Hemoglobin'],
  'hemoglobin-split': ['Examples/Globins/Hemoglobin', 'CopiesSplit'],
  scale: ['Examples/Globins/Rows', 'Scale'],
  'scale-navigators': ['Examples/Globins/Rows', 'ScaleWithNavigators'],
  lazy: ['Examples/Globins/Rows', 'Lazy'],
};

/** The features shown on this data, each a story of its own. */
const Features: readonly (readonly [string, string, string])[] = [
  ['Features/Rows', 'Columns', 'columns'],
  ['Features/Query', 'AlignmentRowsAndColumns', 'query-alignment'],
];

const Editing: readonly (readonly [string, string])[] = [
  ['Family', 'The residue scheme changed (a theme patch) and the logo hidden; Undo and Redo.'],
  ['Folded', 'The folds shown kept (`clades.keep`, one patch), a row pinned above the tree’s order.'],
  ['Hemoglobin', 'The copies split into a block a chain and merged again into ×2.'],
];

const Sources: readonly string[] = [
  'Pfam PF00042 (Globin) seed, Pfam 39.0, from InterPro: CC0 1.0, cut to the 60 rows nearest human hemoglobin (HBA_HUMAN and HBB_HUMAN are not rows of this seed).',
  'UniProtKB, the organism and class of each row: UniProt Consortium, CC BY 4.0.',
  'PDB 4HHB, human deoxyhaemoglobin (Fermi, Perutz, Shaanan, Fourme 1984), wwPDB: CC0 1.0.',
  'The tree (neighbor joining on identity), the heme contacts and the 10,000- and 10,000,000-row variants are computed here.',
];

function item(title: string, exportName: string, label: string, summary: string): HTMLElement {
  const li = document.createElement('li');
  const a = document.createElement('a');
  a.href = `?path=/story/${storyId(title, exportName)}`;
  a.target = '_top';
  a.textContent = label;
  li.append(a, ` ${summary}`);
  return li;
}

function list(items: readonly HTMLElement[]): HTMLElement {
  const ul = document.createElement('ul');
  ul.append(...items);
  return ul;
}

function heading(text: string): HTMLElement {
  const h = document.createElement('h3');
  h.textContent = text;
  return h;
}

/** The Globins views, the features its data shows, its editing views and its sources. */
export const Views: StoryObj = {
  render: () => {
    const root = document.createElement('div');
    root.style.cssText = 'font: 14px/1.5 system-ui, sans-serif; padding: 24px; max-width: 72ch;';
    const h = document.createElement('h2');
    h.textContent = 'Globins';
    const intro = document.createElement('p');
    intro.textContent =
      'The globin family as Pfam’s seed alignment publishes it, its tree, and human hemoglobin (4HHB) with its hemes; generated variants for scale. Each view is one spec over this data.';
    const views = list(GlobinsViews.filter((v) => Where[v.id]).map((v) => item(...Where[v.id], v.title, v.summary)));
    const features = list(
      Features.map(([title, exp, id]) => {
        const v = FeatureViews.find((x) => x.id === id)!;
        return item(title, exp, v.title, v.summary);
      }),
    );
    const editing = list(Editing.map(([exp, summary]) => item('Examples/Globins/Editing', exp, exp, summary)));
    const sources = list(
      Sources.map((t) => {
        const li = document.createElement('li');
        li.textContent = t;
        return li;
      }),
    );
    root.append(
      h,
      intro,
      views,
      heading('Editing'),
      editing,
      heading('Features shown on this data'),
      features,
      heading('Sources'),
      sources,
    );
    return root;
  },
};
