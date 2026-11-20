/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AntibodyViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { storyId } from '../../story-id.js';

const meta: Meta = { title: 'Examples/Antibody/Index' };
export default meta;

/** Where each Antibody view's story is. */
const AntibodyStories: Readonly<Record<string, string>> = {
  schemes: 'Schemes',
  'schemes-design': 'DesignNumberings',
  against: 'DesignAgainstWildType',
  cdrs: 'CdrBookmarks',
  'fc-glycan': 'FcGlycan',
  'fc-glycan-residues': 'GlycanResidues',
  adc: 'Adc',
};

/** The Features shown on this data. */
const Features: readonly (readonly [title: string, summary: string])[] = [
  [
    'Numberings',
    'two declared numberings a chain, read from a source (the wild type) or inline (the design); a ruler in each; a hover’s label; a query by label.',
  ],
  ['Rows', 'the design and the wild type as an alignment of two rows, their indels gaps.'],
  ['Bookmarks', 'a CDR selected by its address, bookmarked from the selection, renamed, recolored, removed.'],
  ['Editing', 'the ruler’s scheme switched; each bookmark edit; every one a step, undone and redone.'],
  [
    'Molecules',
    'an IgG1 Fc (PDB 3AVE) as a complex: its N297 glycan a branched molecule, a chip or its residues (one patch), its links as arcs.',
  ],
];

interface SourceRow {
  readonly file: string;
  readonly entry: string;
  readonly licence: readonly string[];
  readonly attribution: string;
  readonly trimmed: string;
  readonly derived: string;
}

/** The views of the Antibody dataset, its sources, licences and attributions, and what is derived or generated. */
export const Views: StoryObj = {
  render: () => {
    const root = document.createElement('div');
    root.style.cssText =
      'font: var(--msui-text-sm) / 1.5 var(--msui-font-sans); padding: 24px; max-width: 80ch; color: var(--msui-text); background: var(--msui-bg); min-height: 100vh;';
    const el = (tag: 'h2' | 'h3' | 'p', text: string) => {
      const e = document.createElement(tag);
      e.textContent = text;
      return e;
    };
    const views = document.createElement('ul');
    for (const v of AntibodyViews) {
      const li = document.createElement('li');
      const a = document.createElement('a');
      a.href = `?path=/story/${storyId('Examples/Antibody/Views', AntibodyStories[v.id])}`;
      a.target = '_top';
      a.textContent = v.title;
      li.append(a, ` ${v.summary}`);
      views.append(li);
    }
    const features = document.createElement('ul');
    for (const [title, summary] of Features) {
      const li = document.createElement('li');
      li.append(`${title}: ${summary}`);
      features.append(li);
    }
    const sources = document.createElement('ul');
    root.append(
      el('h2', 'Antibody: trastuzumab’s Fab'),
      el(
        'p',
        'PDB 1N8Z, trastuzumab (Herceptin) Fab with HER2: the light chain (A, 214 residues) and the heavy chain (B, 220); HER2 is not taken. PDB 3AVE, a fucosylated human IgG1 Fc with its N297 glycans. Each view is one spec over this data.',
      ),
      views,
      el('h3', 'Features shown on this data'),
      features,
      el('h3', 'Sources, licences and attributions'),
      sources,
      el('h3', 'Generated and derived'),
      el(
        'p',
        'Derived in code: Kabat’s numbering, from SAbDab’s IMGT numbering (counted between the domains’ conserved cysteines and the tryptophan and phenylalanine after them; SAbDab2 serves IMGT’s numbering alone). Definitions in code: the CDRs by Kabat and by IMGT. Generated: the design, from the wild type by an edit list (two substitutions in CDR-H2; a substitution and a two-residue insertion in CDR-H3; a deletion in CDR-L1), numbered by the schemes’ rules; no published design with indels has a permissive licence. Its numbering is not yet compared with ANARCI’s.',
      ),
    );
    // the provenance as recorded, read the way the views read the fixtures
    void fetch('fixtures/antibody/sources.json')
      .then((r) => (r.ok ? r.json() : Promise.reject(new Error(`${r.status}`))))
      .then((rows: SourceRow[]) => {
        for (const r of rows) {
          const li = document.createElement('li');
          // a CC BY source's changes stated (its `trimmed`), as the licence asks
          li.textContent = `${r.file}: ${r.entry}. ${r.licence.join(', ')}; ${r.attribution}. Changed: ${r.trimmed}.${r.derived !== 'none' ? ` Derived: ${r.derived}.` : ''}`;
          sources.append(li);
        }
      })
      .catch((e: unknown) => {
        sources.append(el('p', `The sources could not be read: ${e instanceof Error ? e.message : String(e)}`));
      });
    return root;
  },
};
