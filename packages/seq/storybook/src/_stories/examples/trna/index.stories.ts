/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { TrnaScaleView, TrnaViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { storyId } from '../../story-id.js';

const meta: Meta = { title: 'Examples/tRNA/Index' };
export default meta;

const Exports: Record<string, string> = {
  arcs: 'Arcs',
  split: 'SplitCrossing',
  ring: 'Chords',
  'dot-plot': 'DotPlot',
  'seed-row': 'SeedRow',
  seed: 'Seed',
  'dot-plot-scale': 'DotPlotScale',
};

const Features: readonly (readonly [title: string, summary: string])[] = [
  ['Pairs', 'arcs nested and crossing on one system (1EHZ), chords on a ring, dots on a grid.'],
  [
    'Grid and dot plots',
    'the covariation above the diagonal and the structure’s pairs below, a crosshair, a rectangle.',
  ],
  ['Rows and the logo', 'the seed’s 954 rows under their conservation logo.'],
  [
    'Editing',
    'the tertiary pairs hidden by a filter in the arcs’ `data`; the covariation’s palette changed; each a patch, undone.',
  ],
];

function item(id: string, title: string, summary: string): HTMLElement {
  const li = document.createElement('li');
  const a = document.createElement('a');
  a.href = `?path=/story/${storyId('Examples/tRNA/Views', Exports[id])}`;
  a.target = '_top';
  a.textContent = title;
  li.append(a, ` ${summary}`);
  return li;
}

/** The views of the tRNA dataset, and the features its data shows. */
export const Views: StoryObj = {
  render: () => {
    const root = document.createElement('div');
    root.style.cssText = 'font: 14px/1.5 system-ui, sans-serif; padding: 24px; max-width: 72ch;';
    const h = document.createElement('h2');
    h.textContent = 'tRNA';
    const intro = document.createElement('p');
    intro.textContent =
      'Yeast tRNA-Phe (PDB 1EHZ, chain A: 76 nucleotides, 14 modified residues, 27 base pairs of the entry’s own table) and the Rfam RF00005 seed (954 rows). Each view is one spec over this data; 1EHZ’s place in the seed and the covariation of its positions are computed from the two recorded files by code in the dataset (`derive.ts`) and committed as generated data (`derived-data.ts`, regenerated and checked by a test), so a view builds without reading the seed.';
    const ul = document.createElement('ul');
    ul.append(...[...TrnaViews, TrnaScaleView].map((v) => item(v.id, v.title, v.summary)));
    const h2 = document.createElement('h3');
    h2.textContent = 'Features shown on this data';
    const features = document.createElement('ul');
    for (const [title, summary] of Features) {
      const li = document.createElement('li');
      li.append(`${title}: ${summary}`);
      features.append(li);
    }
    root.append(h, intro, ul, h2, features);
    return root;
  },
};
