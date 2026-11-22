/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { printAddress } from '@molstar/core/repr/tree';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, type SeqSpecBuilder } from '@molstar/seqx-core/lang';
import type { Derived } from './derived.js';
import {
  Cdk2Features,
  constructOffset,
  contactRows,
  crystal,
  CyclinFeatures,
  Data,
  interfaceRows,
  ligandContactRows,
} from './data.js';

// The Complex's crystal, PDB 1FIN: CDK2 and the cyclin A construct, two copies of the complex, an ATP on each CDK2. Its
// chains as blocks of one section, polymers first, then the ATPs as one-position blocks (chips, EC-101); the copies
// merged into one block labeled ×2 (`copies: 'merge'`, the default, EC-65) or split (`each`). Under them the features
// (CDK2's PSTAIRE helix, activation loop and T160; cyclin A's cyclin box folds on the construct), the interface
// residues of both copies, the first copy's secondary structure and the letters; the contacts across the chains and of
// the ATPs with CDK2 as arcs only when asked (`contacts`: optional `pairs` layers, never a default), folded by copy group.

export const title = 'Complex: 1FIN’s chains as blocks';

export interface BlocksOptions {
  readonly arrange?: 'row' | 'stack';
  readonly copies?: 'merge' | 'each';
  readonly layout?: 'rings';
  /** The contacts across the chains as arcs (a `pairs` track). */
  readonly contacts?: boolean;
}

/** The data nodes of 1FIN's views: its sequences, ATP, the features, the interface and the contacts. */
export function crystalData(b: SeqSpecBuilder, d: Derived = Data) {
  const x = crystal(b);
  const atp = b.data.sequence('atp', { alphabet: 'protein', letters: 'X', name: 'ATP' });
  const shift = constructOffset(d);
  const cdk2Features = b.data.ranges('cdk2-features', { target: x.cdk2, rows: Cdk2Features });
  const cyclinFeatures = b.data.ranges('cyclin-features', {
    target: x.cyclin,
    rows: {
      ...CyclinFeatures,
      start: CyclinFeatures.start.map((s) => s - shift),
      end: CyclinFeatures.end.map((e) => e - shift),
    },
  });
  const cdk2Interface = b.data.ranges('cdk2-interface', { target: x.cdk2, rows: interfaceRows('cdk2', d) });
  const cyclinInterface = b.data.ranges('cyclin-interface', { target: x.cyclin, rows: interfaceRows('cyclin', d) });
  const contacts = b.data.pairs('contacts', { a: x.cdk2, b: x.cyclin, rows: contactRows({}, d) });
  // each ATP's contacts with CDK2's residues, both copies' in one node (the ATP block is one position, the copies' merged)
  const atps = d.ligands.filter((l) => l.entry === '1FIN');
  const atpContacts = b.data.pairs('atp-contacts', {
    a: atp,
    b: x.cdk2,
    rows: ligandContactRows({ ...atps[0], contacts: atps.flatMap((l) => l.contacts) }),
  });
  return { ...x, atp, cdk2Features, cyclinFeatures, cdk2Interface, cyclinInterface, contacts, atpContacts };
}

/** A layer's data: a ranges node's rows of one chain (`(item :chain A)`). */
export const ofChain = (node: string, chain: string) => `${node}${printAddress([{ kind: 'item', params: { chain } }])}`;

/** The complex's chains and ATPs as blocks, its features, interface and letters under them. */
export function build(options: BlocksOptions = {}, d: Derived = Data): Node {
  const b = createSeqSpec({ title });
  const x = crystalData(b, d);
  b.views.view('main', (v) => {
    const s = v.section('complex', {
      ...(options.arrange ? { arrange: options.arrange } : {}),
      ...(options.copies ? { copies: options.copies } : {}),
      ...(options.layout ? { layout: options.layout } : {}),
    });
    for (const c of d.chains)
      s.block(c.id.toLowerCase(), {
        of: c.system === 'cdk2' ? x.cdk2 : x.cyclin,
        label: `${c.system === 'cdk2' ? 'CDK2' : 'cyclin A'} (${c.id})`,
      });
    for (const l of d.ligands.filter((l) => l.entry === '1FIN'))
      s.block(`atp-${l.chain.toLowerCase()}`, { of: x.atp, label: 'ATP' });
    s.track('features', { label: 'Features' }, (t) => {
      t.layer('cdk2', { type: 'ranges', data: x.cdk2Features });
      t.layer('cyclin', { type: 'ranges', data: x.cyclinFeatures });
    });
    s.track('interface', { label: 'Interface' }, (t) => {
      t.layer('cdk2', { type: 'ranges', data: x.cdk2Interface });
      t.layer('cyclin', { type: 'ranges', data: x.cyclinInterface });
    });
    if (options.contacts)
      s.track('contacts', { label: 'Contacts (≤ 4 Å)' }, (t) => {
        t.layer('arcs', { type: 'pairs', data: x.contacts });
        t.layer('atp', { type: 'pairs', data: x.atpContacts });
      });
    s.track('secondary', { label: 'Secondary (A, B)' }, (t) => {
      t.layer('cdk2', { type: 'secondary', data: ofChain(x.cdk2Secondary, 'A') });
      t.layer('cyclin', { type: 'secondary', data: ofChain(x.cyclinSecondary, 'B') });
    });
    s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** The copies split: a block a chain (`copies: 'each'`). */
export const split = () => build({ copies: 'each' });

/** The chains stacked: a block a row, aligned at their first residue (`arrange: 'stack'`). */
export const stacked = () => build({ arrange: 'stack' });

/** The chains stacked with the copies split: a row a chain, the two ATPs as chips on one row of their own. */
export const stackedSplit = () => build({ arrange: 'stack', copies: 'each' });

/** The chains on a ring: a block an arc. */
export const ring = () => build({ layout: 'rings' });

/** The contacts across the chains, and each ATP with CDK2, as arcs between the merged blocks: one arc a contact of either copy. */
export const contacts = () => build({ contacts: true });

/** The contacts with the copies split: each copy's contacts between its own blocks. */
export const contactsSplit = () => build({ contacts: true, copies: 'each' });

/** The copies' scale variant (performance.md, Scale): `n` copies of CDK2 as blocks, merged into one (`×n`) or each. */
export function copiesScale(n = 200, copies: 'merge' | 'each' = 'merge'): Node {
  const b = createSeqSpec({ title: `Complex: ${n} copies of CDK2` });
  const x = crystal(b);
  const features = b.data.ranges('cdk2-features', { target: x.cdk2, rows: Cdk2Features });
  b.views.view('main', (v) => {
    const s = v.section('copies', { copies });
    for (let i = 0; i < n; i++) s.block(`c${i + 1}`, { of: x.cdk2, label: `CDK2 ${i + 1}` });
    s.track('features', { label: 'Features' }).layer('cdk2', { type: 'ranges', data: features });
    s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}
