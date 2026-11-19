/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch, type Patch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { addRefPatch, createSeqSpec, type RefInput, SeqLanguage } from '@molstar/seqx-core/lang';
import { HemoglobinFile } from './data.js';
import { HemeContacts } from './hemes.js';

// Hemoglobin (examples.md, Globins; EC-64, EC-65): PDB 4HHB, α2β2, as blocks of one section: the α chains A and C
// and the β chains B and D (two molecules, four units: a copy group each), then the four hemes, each a one-position
// sequence drawn as a chip; the chains' helices A to H from the entry; each heme's contacts (heavy atoms within 4 Å,
// derived) as arcs to its chain and as the "heme pocket" ranges, bookmarked. The copies merged into one block labeled
// ×2 (`copies: 'merge'`, the default) or split (`each`), a patch each way. The heme per heavy atom and the units' link
// to their rows in the family wait: a heme is one position this slice, and the language has no alignment map.

export const title = 'Hemoglobin (4HHB)';

/**
 * The chains' lengths (anchors: 141 and 146 residues, the heme one position): each block names its range, as a section
 * of blocks reads its systems' lengths when it resolves, before a source has said them.
 */
export const HemoglobinLengths = { hba: 141, hbb: 146, hem: 1 } as const;

/** The section's address and its path (what a patch names). */
export const HemoglobinSection = '/views/main/hb';
const SectionIds = ['views', 'main', 'hb'];

/** The pocket's ranges on a chain: its contacting residues, 0-based, one range a residue. */
const pocket = (residues: readonly number[]) => ({
  start: residues.map((r) => r - 1),
  end: residues.map((r) => r),
  id: residues.map((r) => `r${r}`),
});

/** Contacts from the heme (its one position) to its chain's residues. */
const contacts = (residues: readonly number[]) => ({
  a_start: residues.map(() => 0),
  a_end: residues.map(() => 1),
  b_start: residues.map((r) => r - 1),
  b_end: residues.map((r) => r),
  type: residues.map(() => 'contact'),
});

/** The bookmarks the builder writes: each molecule's heme pocket. */
export const PocketRefs: readonly RefInput[] = [
  { id: 'pocket-a', to: '/data/hba-pocket', label: 'Heme pocket, α' },
  { id: 'pocket-b', to: '/data/hbb-pocket', label: 'Heme pocket, β' },
];

export interface HemoglobinOptions {
  readonly copies?: 'merge' | 'each';
}

/** 4HHB's chains and hemes as blocks, their helices, pockets and heme contacts under them. */
export function build(options: HemoglobinOptions = {}): Node {
  const b = createSeqSpec({ title });
  const file = b.data.source('4hhb', {
    format: 'seq',
    uri: HemoglobinFile,
    attribution: 'wwPDB, PDB 4HHB (Fermi, Perutz, Shaanan, Fourme 1984)',
    licence: 'CC0-1.0',
  });
  const hba = b.data.sequence('hba', { alphabet: 'protein', source: file, name: 'hba' });
  const hbb = b.data.sequence('hbb', { alphabet: 'protein', source: file, name: 'hbb' });
  const hem = b.data.sequence('hem', { alphabet: 'protein', source: file, name: 'hem' });
  const helicesA = b.data.ranges('hba-helices', { target: hba, source: file });
  const helicesB = b.data.ranges('hbb-helices', { target: hbb, source: file });
  const pocketA = b.data.ranges('hba-pocket', { target: hba, rows: pocket(HemeContacts.hba) });
  const pocketB = b.data.ranges('hbb-pocket', { target: hbb, rows: pocket(HemeContacts.hbb) });
  const contactsA = b.data.pairs('hba-contacts', { a: hem, b: hba, rows: contacts(HemeContacts.hba) });
  const contactsB = b.data.pairs('hbb-contacts', { a: hem, b: hbb, rows: contacts(HemeContacts.hbb) });
  b.views.view('main', (v) => {
    const s = v.section('hb', options.copies ? { copies: options.copies } : {});
    const L = HemoglobinLengths;
    s.block('a', { of: hba, label: 'α (A)', range: [0, L.hba] });
    s.block('b', { of: hbb, label: 'β (B)', range: [0, L.hbb] });
    s.block('c', { of: hba, label: 'α (C)', range: [0, L.hba] });
    s.block('d', { of: hbb, label: 'β (D)', range: [0, L.hbb] });
    for (const c of ['a', 'b', 'c', 'd']) s.block(`hem-${c}`, { of: hem, label: 'heme', range: [0, L.hem] });
    s.track('contacts', { label: 'Heme contacts' }, (t) => {
      t.layer('a', { type: 'pairs', data: contactsA });
      t.layer('b', { type: 'pairs', data: contactsB });
    });
    s.track('helices', { label: 'Helices' }, (t) => {
      t.layer('a', { type: 'ranges', data: helicesA });
      t.layer('b', { type: 'ranges', data: helicesB });
    });
    s.track('pocket', { label: 'Heme pocket' }, (t) => {
      t.layer('a', { type: 'ranges', data: pocketA });
      t.layer('b', { type: 'ranges', data: pocketB });
    });
    s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
  });
  let root = b.build();
  for (const ref of PocketRefs) {
    const r = applyPatch(root, addRefPatch(root, ref), { language: SeqLanguage });
    if (!r.ok) throw new Error(`The bookmark ${ref.id} did not apply: ${r.issues.map((i) => i.code).join(', ')}`);
    root = r.root;
  }
  return root;
}

/** The copies split: a block a chain and a heme each. */
export const split = () => build({ copies: 'each' });

/** The editing view's toggle: the copies merged or split, one patch. */
export function copiesPatch(copies: 'merge' | 'each'): Patch {
  return { ops: [{ op: 'set', at: SectionIds, key: 'copies', value: copies }] };
}
