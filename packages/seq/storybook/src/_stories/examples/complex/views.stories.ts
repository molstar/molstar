/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ComplexViews, type ExampleView } from '@molstar/seqx-examples';
import { copiesScale } from '@molstar/seqx-examples/complex/blocks';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';
import { LinksHoverProvider } from '@molstar/seqx-viewer/behaviors';
import { hideLigandsPast, setCopies, setPlddtBins, toggleAtpForm } from './compose.js';
import { complexEditStory } from './edit-root.js';

const view = (id: string) => ComplexViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/Complex/Views' };
export default meta;

/**
 * Both chains of the prediction with AlphaFold DB's per-residue confidence in its four bins, UniProt's secondary
 * structure and 1FIN's mapped onto UniProt's numbering, and what 1FIN holds: cyclin A2's first 172 residues gray (not
 * in the crystal). The editing view: other bins for the confidence (one theme patch), undone.
 */
export const Confidence: StoryObj = {
  render: () =>
    complexEditStory({
      view: view('confidence'),
      edits: [
        { label: 'Bins at 60, 80', run: (s) => setPlddtBins(s, [0, 60, 80, 90, 100]) },
        { label: 'AlphaFold’s bins', run: (s) => setPlddtBins(s, [0, 50, 70, 90, 100]) },
      ],
    }),
};

/**
 * CDK2 and its ligands as blocks after it: 1FIN's ATP, then the ligands of six CC0 entries of CDK2, best resolution
 * first, each a chip with arcs to its contacting residues. Hover an arc: the ligand and its residue light; run a
 * ligand's contacts from the bookmarks: its residues highlighted. The editing view: hide the ligands of entries past
 * 2 Å (a filter in each layer's data, one patch), undone.
 */
export const Ligands: StoryObj = {
  render: () =>
    complexEditStory({
      view: view('ligands'),
      edits: [
        { label: 'Hide ligands past 2 Å', run: (s) => hideLigandsPast(s, 2) },
        { label: 'Show every ligand', run: (s) => hideLigandsPast(s, Infinity) },
      ],
    }),
};

/** ATP as its 31 heavy atoms, a block of its own: each atom's contacts an arc, its name under it. */
export const AtpPerAtom: StoryObj = { render: () => exampleStory(view('ligands-atoms')) };

/** 1FIN's first copy as text: CDK2, then cyclin A, then ATP as a chip flowing with the letters like a word. */
export const Text: StoryObj = { render: () => exampleStory(view('text')) };

/** The copies' scale: 200 copies of CDK2 as blocks, merged into one labeled ×200, split by one patch (EC-65). */
const scale: ExampleView = {
  dataset: 'complex',
  id: 'copies-scale',
  title: '200 copies',
  summary:
    'The copies’ scale variant (generated): 200 copies of CDK2, merged into one block labeled ×200; split them (one patch) and every copy is a block of its own.',
  build: () => copiesScale(),
};

export const Copies: StoryObj = {
  render: () =>
    complexEditStory({
      view: scale,
      edits: [
        { label: 'Split the copies', run: (s) => setCopies(s, 'copies', 'each') },
        { label: 'Merge the copies', run: (s) => setCopies(s, 'copies', 'merge') },
      ],
    }),
};

/** The host's choice the stories show on (off in the default viewer): a ligand lights its contacts, a residue its ligand. */
const LinksHoverOn = [{ provider: LinksHoverProvider, params: undefined }] as const;

/**
 * PDB 1FIN as a complex of its units: CDK2 ×2 and cyclin A ×2 with their author numbers (the construct from 173), ATP
 * ×2 a chip colored on its own, then the six CC0 entries' ligands (standing for AlphaFill's transplants) with their
 * resolution as a score; the contacts per pair of units as arcs. The transplants and their contacts come from those
 * other entries (computed against each entry's CDK2, drawn on 1FIN's chain A): a transplant's card and its arcs name
 * its entry. Hover ATP's chip: the CDK2 residues it contacts light (`links.hover`, on here); hover a contact: both
 * ends. Editing: ATP as a chip or by its heavy atoms (one patch).
 */
export const AsAComplex: StoryObj = {
  render: () =>
    complexEditStory({
      view: view('units'),
      edits: [{ label: 'ATP: chip/atoms', run: (s) => toggleAtpForm(s) }],
      behaviors: LinksHoverOn,
    }),
};

/** The complex with the first ATP by its heavy atoms, each labeled by its element. */
export const AsAComplexAtpAtoms: StoryObj = {
  render: () => exampleStory(view('units-atoms'), { behaviors: LinksHoverOn }),
};
