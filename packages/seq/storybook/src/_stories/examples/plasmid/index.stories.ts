/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PlasmidViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Figures, idOf, Tries, type Where } from './where.js';

const meta: Meta = { title: 'Examples/Plasmid/Index' };
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

/** The views of the Plasmid dataset: its figures, then the editing views and the Features that act on its data. */
export const Views: StoryObj = {
  render: () => {
    const root = document.createElement('div');
    root.style.cssText = 'font: 14px/1.5 system-ui, sans-serif; padding: 24px; max-width: 72ch;';
    const h = document.createElement('h2');
    h.textContent = 'Plasmid';
    const intro = document.createElement('p');
    intro.textContent =
      'pUC19c, the cloning vector (INSDC L09137.2, 2,686 bp, circular). The record has no feature annotation, so its features are found in the sequence in code (bla and lacZα as the open reading frames that encode their proteins, the cloning site from EcoRI to HindIII, the pMB1 origin by its signature) and its restriction sites computed on the cyclic sequence. The rotated variant moves the origin inside bla, so bla and a generated primer site cross it.';
    const licence = document.createElement('p');
    licence.textContent =
      'Source: INSDC L09137.2 (ENA, NCBI); Yanisch-Perron C, Vieira J, Messing J, Gene 33:103-119 (1985). NCBI places no restriction on the record (LicenseRef-NCBI-public-data); the features, sites and variants are generated in code (LicenseRef-seqx-generated).';
    const figures = list(PlasmidViews.filter((v) => Figures[v.id]).map((v) => item(Figures[v.id], v.title, v.summary)));
    const h2 = document.createElement('h3');
    h2.textContent = 'Editing, and Features shown on this data';
    const features = list(Tries.map((t) => item(t.where, t.title, t.summary)));
    root.append(h, intro, licence, figures, h2, features);
    return root;
  },
};
