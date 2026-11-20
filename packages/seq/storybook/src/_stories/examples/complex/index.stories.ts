/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ComplexViews } from '@molstar/seqx-examples';
import { AlphaFillPolicy } from '@molstar/seqx-examples/complex/data';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { storyId } from '../../story-id.js';

const meta: Meta = { title: 'Examples/Complex/Index' };
export default meta;

/** Where each Complex view's story is. */
const ComplexStories: Readonly<Record<string, readonly [title: string, name: string]>> = {
  blocks: ['Examples/Complex/Blocks', 'Row'],
  'blocks-split': ['Examples/Complex/Blocks', 'Split'],
  'blocks-stacked': ['Examples/Complex/Blocks', 'Stacked'],
  'blocks-stacked-split': ['Examples/Complex/Blocks', 'StackedSplit'],
  'blocks-ring': ['Examples/Complex/Blocks', 'OnARing'],
  contacts: ['Examples/Complex/PAE and contacts', 'Contacts'],
  pae: ['Examples/Complex/PAE and contacts', 'Pae'],
  'complex-grid': ['Examples/Complex/PAE and contacts', 'OverTheComplex'],
  confidence: ['Examples/Complex/Views', 'Confidence'],
  ligands: ['Examples/Complex/Views', 'Ligands'],
  'ligands-atoms': ['Examples/Complex/Views', 'AtpPerAtom'],
  text: ['Examples/Complex/Views', 'Text'],
  units: ['Examples/Complex/Views', 'AsAComplex'],
  'units-atoms': ['Examples/Complex/Views', 'AsAComplexAtpAtoms'],
};

/** The Features shown on this data, each with what to try. */
const Features: readonly (readonly [title: string, name: string, label: string, summary: string])[] = [
  ['Examples/Complex/Blocks', 'Editing', 'Editing', 'Bookmark the interface; merge and split the copies; undo, redo.'],
  [
    'Examples/Complex/Views',
    'Copies',
    'Copies at scale',
    '200 copies of CDK2 merged into one block, split by a patch.',
  ],
  ['Examples/Complex/Views', 'Confidence', 'Themes', 'The confidence’s bins changed by one theme patch.'],
  [
    'Examples/Complex/Views',
    'Ligands',
    'Filters',
    'The ligands of entries past 2 Å hidden by a filter in each layer’s data.',
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

function link(title: string, name: string, text: string): HTMLAnchorElement {
  const a = document.createElement('a');
  a.href = `?path=/story/${storyId(title, name)}`;
  a.target = '_top';
  a.textContent = text;
  return a;
}

/** The views of the Complex dataset, its sources, licences and attributions, and what is derived or generated. */
export const Views: StoryObj = {
  render: () => {
    const root = document.createElement('div');
    root.style.cssText =
      'font: var(--msui-text-sm) / 1.5 var(--msui-font-sans); padding: 24px; max-width: 80ch; color: var(--msui-text); background: var(--msui-bg); min-height: 100vh;';
    const h = (level: 'h2' | 'h3', text: string) => {
      const el = document.createElement(level);
      el.textContent = text;
      return el;
    };
    const p = (text: string) => {
      const el = document.createElement('p');
      el.textContent = text;
      return el;
    };
    const views = document.createElement('ul');
    for (const v of ComplexViews) {
      const where = ComplexStories[v.id];
      const li = document.createElement('li');
      li.append(where ? link(where[0], where[1], v.title) : v.title, ` ${v.summary}`);
      views.append(li);
    }
    const features = document.createElement('ul');
    for (const [title, name, label, summary] of Features) {
      const li = document.createElement('li');
      li.append(link(title, name, label), ` ${summary}`);
      features.append(li);
    }
    const sources = document.createElement('ul');
    root.append(
      h('h2', 'Complex: CDK2 with cyclin A'),
      p(
        'The crystal PDB 1FIN (two copies of the complex, an ATP on each CDK2) and the ligands of six CC0 entries of CDK2; the prediction of each chain from AlphaFold DB, their UniProt entries. Each view is one spec over this data.',
      ),
      views,
      h('h3', 'Features shown on this data'),
      features,
      h('h3', 'Sources, licences and attributions'),
      sources,
      h('h3', 'Generated and derived'),
      p(
        `Derived at recording, from the coordinates: the contacts within 4 Å between heavy atoms across each copy’s chains, and each ligand’s contacting residues. Generated in code: the PAE over the complex (synthetic: no permissive prediction of the complex was found; low across 1FIN’s interface, from a seed), and the 200 copies of the scale view. Curated in code: the PSTAIRE helix, the activation loop and the cyclin box folds (UniProt has no feature of them). AlphaFill’s transplants are not taken: its usage policy (${AlphaFillPolicy.url}, checked ${AlphaFillPolicy.checked}) states no CC BY 4.0 or looser licence.`,
      ),
    );
    // the provenance as recorded, read the way the views read the fixtures
    void fetch('fixtures/complex/sources.json')
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
        sources.append(p(`The sources could not be read: ${e instanceof Error ? e.message : String(e)}`));
      });
    return root;
  },
};
