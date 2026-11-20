/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch, type Patch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { addRefPatch, createSeqSpec, type RefInput, SeqLanguage } from '@molstar/seqx-core/lang';
import { HemoglobinFile } from './data.js';
import { HemeContacts, HemeUnitContacts, HemoglobinChains } from './hemes.js';

// Hemoglobin (examples.md, Globins; EC-64, EC-65): PDB 4HHB, α2β2, a complex of eight units: the α chains A and C
// and the β chains B and D (two sequences, a copy group each), and the four hemes, each a unit of a one-position
// sequence drawn as a chip; blocks of one section in that order. The chains' helices A to H from the entry; the hemes'
// contacts (heavy atoms within 4 Å, derived per pair of units) as the complex's pairs between units, so each heme's
// arcs reach the chain it sits in; each molecule's heme pocket as ranges, bookmarked. The copies merged into one block
// labeled ×2 (`copies: 'merge'`, the default) or split (`each`), a patch each way. The heme per heavy atom and the
// units' link to their rows in the family wait: a heme is one position this slice, and the language has no alignment
// map.

export const title = 'Hemoglobin (4HHB)';

/** The chains' lengths (anchors: 141 and 146 residues, the heme one position): the file's, which the blocks read. */
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

/** The hemes' contacts as the complex's pairs: a row a residue, from the heme's unit (its one position) to its chain's. */
function contacts() {
  const rows = HemeUnitContacts.flatMap((c) => c.residues.map((r) => ({ heme: hemeUnit(c.heme), chain: c.chain, r })));
  return {
    a_unit: rows.map((x) => x.heme),
    a_start: rows.map(() => 0),
    a_end: rows.map(() => 1),
    b_unit: rows.map((x) => x.chain),
    b_start: rows.map((x) => x.r - 1),
    b_end: rows.map((x) => x.r),
    type: rows.map(() => 'contact'),
  };
}

/** A heme's unit id by its chain: `hem-a` … */
const hemeUnit = (chain: string) => `hem-${chain.toLowerCase()}`;

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
  const sequences = { hba, hbb };
  const units: Record<string, string> = {};
  b.data.complex('hb', { name: 'Hemoglobin (4HHB)' }, (c) => {
    for (const ch of HemoglobinChains)
      units[ch.chain] = c.unit(ch.chain, { of: sequences[ch.of], chain: ch.chain, auth_chain: ch.chain });
    for (const ch of HemoglobinChains)
      units[hemeUnit(ch.chain)] = c.unit(hemeUnit(ch.chain), { of: hem, chain: ch.heme, auth_chain: ch.chain });
    c.links('contacts', { rows: contacts() });
  });
  b.views.view('main', (v) => {
    const s = v.section('hb', options.copies ? { copies: options.copies } : {});
    // each block its whole chain: the lengths the file gives once it is open
    for (const ch of HemoglobinChains)
      s.block(ch.chain.toLowerCase(), {
        of: units[ch.chain] as never,
        label: `${ch.of === 'hba' ? 'α' : 'β'} (${ch.chain})`,
      });
    for (const ch of HemoglobinChains)
      s.block(hemeUnit(ch.chain), { of: units[hemeUnit(ch.chain)] as never, label: 'heme' });
    s.track('contacts', { label: 'Heme contacts' }, (t) =>
      t.layer('contacts', { type: 'pairs', data: '/data/hb/contacts' as never }),
    );
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
