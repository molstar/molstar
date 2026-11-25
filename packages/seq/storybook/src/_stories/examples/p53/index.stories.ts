/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { FeatureStories, FeatureTries, Figures, idOf, type Where } from './where.js';

const meta: Meta = { title: 'Examples/p53/Index' };
export default meta;

function item(href: Where, title: string, summary: string): HTMLElement {
  const li = document.createElement('li');
  const a = document.createElement('a');
  a.href = `?path=/story/${idOf(href)}`;
  a.target = '_top';
  a.textContent = title;
  li.append(a, ` ${summary}`);
  return li;
}

function list(items: readonly HTMLElement[]): HTMLElement {
  const ul = document.createElement('ul');
  ul.append(...items);
  return ul;
}

/** The views of the p53 dataset: its figures, then the features its data shows (each a story of its own). */
export const Views: StoryObj = {
  render: () => {
    const root = document.createElement('div');
    root.style.cssText = 'font: 14px/1.5 system-ui, sans-serif; padding: 24px; max-width: 72ch;';
    const h = document.createElement('h2');
    h.textContent = 'p53';
    const intro = document.createElement('p');
    intro.textContent =
      'Human TP53: its protein with domains and conservation, and a toy gene it is placed on, which encodes it but for one variant (R175H). Each view is one spec over this data.';
    const figures = list(P53Views.filter((v) => Figures[v.id]).map((v) => item(Figures[v.id], v.title, v.summary)));
    const h2 = document.createElement('h3');
    h2.textContent = 'Features shown on this data';
    const note = document.createElement('p');
    note.textContent = 'How each works is a Features story, on whichever dataset shows it best; these are on p53.';
    const features = list([
      ...P53Views.filter((v) => FeatureStories[v.id]).map((v) => item(FeatureStories[v.id], v.title, v.summary)),
      ...FeatureTries.map((t) => item(t.where, t.title, t.summary)),
    ]);
    root.append(h, intro, figures, h2, note, features);
    return root;
  },
};
