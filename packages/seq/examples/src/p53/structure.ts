/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, identifiers } from '@molstar/seqx-core/lang';
import { addDomains, P53 } from './data.js';

// p53's structure (examples.md, p53; structure-syncing.md, Datasets): PDB 2OCJ, four copies of the DNA-binding domain
// (P04637 94–312) with no DNA, a zinc in each. The entity's letters and chain A's secondary structure from its `seq`
// file; a complex of its units: the chains A–D (modeled 96–289, the ends unobserved, gray; their `auth` numbering
// UniProt's) and the zincs E–H (`ZN`, author numbers 501–504), each zinc's four metal links (C176, H179, C238, C242)
// the complex's. P04637 onto the entity by a `structure` map (SIFTS's one segment, ranked first) and by a `segments`
// map of the same pair (the global alignment of the two, the stopgap of an `alignment` map, ranked lower): the
// section's paths (P04637's domains drawn on the chains) take either by `via`. The units, their observed runs, their author numbers and the zinc sites are the entry's, held
// here and checked against the committed chain table (`2OCJ-chains.json`, the fake structure viewer's own) by a test.

export const title = 'p53 and PDB 2OCJ';

/** The entity's `seq` file (`scripts/record.ts p53-structure`) and the fake structure viewer's chain table. */
export const StructureFile = 'fixtures/p53/2OCJ.seq';
export const ChainTableFile = 'fixtures/p53/2OCJ-chains.json';

/** The entity on P04637: positions 93 to 312 (0-based, 219 residues). */
export const EntityStart = 93;
export const EntityLength = 219;
/** Each chain's modeled residues on the entity, 0-based, end exclusive (`label_seq_id` 3 to 196: 96 to 289). */
export const Observed: readonly [number, number] = [2, 196];
/** The chains (`label_asym_id` = `auth_asym_id`) and each one's zinc: its `label_asym_id` and author number. */
export const Chains = [
  { chain: 'A', zinc: 'E', zincAuth: '501' },
  { chain: 'B', zinc: 'F', zincAuth: '502' },
  { chain: 'C', zinc: 'G', zincAuth: '503' },
  { chain: 'D', zinc: 'H', zincAuth: '504' },
] as const;
/** A zinc's ligands (`struct_conn`, `metalc`): the residue's author number and atom. */
export const ZincSites = [
  { auth: 176, atom: 'SG' },
  { auth: 179, atom: 'ND1' },
  { auth: 238, atom: 'SG' },
  { auth: 242, atom: 'SG' },
] as const;
/** P04637 on the entity by their global alignment (`alignedSegments`, checked by a test): one run, no gap. */
export const AlignedSegments = { from: [EntityStart], to: [0], length: [EntityLength] } as const;

/** The entity position of an author number: `auth` is P04637's number. */
export const positionOf = (auth: number) => auth - EntityStart - 1;

export interface StructureOptions {
  readonly copies?: 'merge' | 'each';
  /** The map a section's paths take (`via`): SIFTS's or the alignment's; absent, the ranked one (SIFTS). */
  readonly via?: 'sifts' | 'aln';
}

/** The addresses a host's adapter and the stories name. */
export const StructureAddresses = {
  protein: '/data/p53',
  entity: '/data/2ocj-1',
  complex: '/data/2ocj',
  sifts: '/data/sifts',
  aln: '/data/aln',
  section: '/views/main/structure',
} as const;

/**
 * P04637 and 2OCJ's units: P04637's domains and letters, then a section of the complex (copies merged by default) with
 * the zincs' links, P04637's domains on the chains, the entity's secondary structure and letters.
 */
export function build(options: StructureOptions = {}): Node {
  const b = createSeqSpec({ title });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53, name: 'P04637' });
  const domains = addDomains(b, p53);
  const file = b.data.source('2ocj-file', {
    format: 'seq',
    uri: StructureFile,
    attribution: 'wwPDB, PDB 2OCJ (Wang, Rosengarth, Luecke 2007)',
    licence: 'CC0-1.0',
  });
  const entity = b.data.sequence('2ocj-1', { alphabet: 'protein', source: file, name: '2OCJ-1' });
  const secondary = b.data.ranges('2ocj-secondary', { target: entity, source: file });
  const zn = b.data.molecule('zn', {
    class: 'ion',
    name: 'zinc',
    residues: ['ZN'],
    atoms: { residue: [0], name: ['ZN'], element: ['ZN'] },
    identifiers: identifiers([['ccd', 'ZN']]),
  });
  const labels = Array.from({ length: EntityLength }, (_, i) =>
    i < Observed[0] || i >= Observed[1] ? null : String(EntityStart + 1 + i),
  );
  b.data.complex('2ocj', { name: '2OCJ' }, (c) => {
    for (const { chain } of Chains)
      c.unit(
        chain,
        { of: entity, chain, auth_chain: chain, observed: { start: [Observed[0]], end: [Observed[1]] } },
        (u) => u.numbering('auth', { method: 'labels', labels }),
      );
    for (const { chain, zinc, zincAuth } of Chains)
      c.unit(zinc, { of: zn, chain: zinc, auth_chain: chain, label: 'Zn' }, (u) =>
        u.numbering('auth', { method: 'labels', labels: [zincAuth] }),
      );
    const rows = Chains.flatMap(({ chain, zinc }) => ZincSites.map((s) => ({ chain, zinc, ...s })));
    c.links('links', {
      rows: {
        a_unit: rows.map((r) => r.chain),
        a_start: rows.map((r) => positionOf(r.auth)),
        a_end: rows.map((r) => positionOf(r.auth) + 1),
        b_unit: rows.map((r) => r.zinc),
        b_start: rows.map(() => 0),
        b_end: rows.map(() => 1),
        type: rows.map(() => 'metal'),
        a_atom: rows.map((r) => r.atom),
        b_atom: rows.map(() => 'ZN'),
      },
    });
  });
  b.data.map('sifts', {
    from: p53,
    to: entity,
    method: 'structure',
    segments: { from: [EntityStart], to: [0], length: [EntityLength], provenance: ['sifts'] },
    priority: 1,
  });
  b.data.map('aln', {
    from: p53,
    to: entity,
    method: 'segments',
    segments: { from: [...AlignedSegments.from], to: [...AlignedSegments.to], length: [...AlignedSegments.length] },
  });
  b.views.view('main', (v) => {
    const p = v.section('protein', { of: p53 });
    p.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    p.track('domains', { label: 'Domains' }).layer('domains', { type: 'ranges', data: domains });
    p.track('letters', { label: 'P04637' }).layer('letters', { type: 'letters' });
    const s = v.section('structure', {
      of: StructureAddresses.complex as never,
      ...(options.copies ? { copies: options.copies } : {}),
      ...(options.via ? { via: StructureAddresses[options.via] as never } : {}),
    });
    s.track('links', { label: 'Zinc' }).layer('links', {
      type: 'pairs',
      data: `${StructureAddresses.complex}/links` as never,
    });
    // P04637's domains on the chains: through the map the section takes (SIFTS's, or the one `via` names)
    s.track('domains', { label: 'Domains' }).layer('domains', { type: 'ranges', data: domains });
    s.track('secondary', { label: 'Secondary structure' }).layer('secondary', { type: 'ranges', data: secondary });
    s.track('letters', { label: '2OCJ' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** The copies split: a block a chain, a chip a zinc. */
export const split = () => build({ copies: 'each' });

/** The section's paths through the alignment's map, not SIFTS's. */
export const viaAlignment = () => build({ via: 'aln' });
