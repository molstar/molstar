/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A complex's unit facts (_spec/README.md "Structures", Facts): the author numbering both ways (EC-1 to EC-8), copies
// (EC-65), observed runs (EC-118), residue rows; and the round trip a description makes through `structureNodes` and
// `structureFacts`: every chain's ids, its molecule and kind, each residue's author number, insertion code and observed
// flag as described, on the zinc finger and 1EHZ (quick) through the mmCIF runtime's own reading of the file.

import type { Node } from '@molstar/core/repr/tree';
import { type NodeEntry, nodeKey } from '@molstar/seqx-core/data';
import { type StructureDescription, structureNodes } from '@molstar/seqx-core/io/structure';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { buildOverlay, nodeAddress, resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { Source, zincFinger } from '../../io/_test/fixtures.js';
import { mmcifDescription } from '../../runtimes/mmcif-description.js';
import { type EntryId, frameOf, openEntry } from '../../runtimes/_test/mmcif-entries.js';
import { structureFacts, type UnitFacts } from '../facts.js';
import type { StructureModel } from '../model.js';
import { chainId, Complex, Ent, EntLength, EntStart, Observed, p53Spec } from './p53.js';

describe('structureFacts', () => {
  const model = resolveSpec(p53Spec({ entityAuth: true }));
  const facts = structureFacts(model, Complex)!;
  const a = facts.units[0];

  it('per unit its ids, molecule, kind, range and copy group (EC-65)', () => {
    expect(facts.units.map((u) => [u.id, u.chain, u.authChain, u.kind, u.class])).toEqual([
      ...[0, 1, 2, 3].map((k) => [chainId(k), chainId(k), chainId(k), 'polymer', 'polymer']),
      ...[0, 1, 2, 3].map((k) => [`Z${chainId(k)}`, `Z${chainId(k)}`, chainId(k), 'other', 'ion']),
    ]);
    expect([a.molecule, a.start, a.length, Array.from(a.observed)]).toEqual([Ent, 0, EntLength, [...Observed]]);
    expect(new Set(facts.units.slice(0, 4).map((u) => u.copy)).size).toBe(1);
    expect(facts.units[4].copy).not.toBe(a.copy);
    expect(Array.from(facts.units[4].observed)).toEqual([0, 1]);
    expect(structureFacts(model, '/data/none')).toBeUndefined();
  });

  it('the unit’s own author numbering beats its molecule’s of the same id (EC-8); both ways', () => {
    expect(a.auth!.label(Observed[0])).toBe(String(EntStart + 1 + Observed[0]));
    expect(a.auth!.position(String(EntStart + 1 + 50))).toBe(50);
    // an insertion code is part of the label (EC-5); an unmodeled residue has none (EC-1 to EC-8 at the edge)
    expect([a.auth!.label(10), a.auth!.position(`${EntStart + 11}A`), a.auth!.position(`${EntStart + 11}a`)]).toEqual([
      `${EntStart + 11}A`,
      10,
      10,
    ]);
    expect([a.auth!.label(0), a.auth!.label(-1), a.auth!.label(EntLength), a.auth!.position('99999')]).toEqual([
      undefined,
      undefined,
      undefined,
      -1,
    ]);
  });

  it('a unit with no numbering of its own reads its molecule’s, shifted onto it (EC-8)', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'ACDEFGHIK' }, (q) =>
      q.numbering('auth', { method: 'labels', labels: ['10', '11', '12', '12A', '13', '14', '15', '16', '17'] }),
    );
    b.data.complex('c', {}, (c) => c.unit('A', { of: s, range: [2, 6] }));
    const u = structureFacts(resolveSpec(b.build()), '/data/c')!.units[0];
    expect([u.start, u.length, u.auth!.label(0), u.auth!.label(1), u.auth!.position('13')]).toEqual([
      2,
      4,
      '12',
      '12A',
      2,
    ]);
    expect(u.auth!.position('10')).toBe(-1);
  });

  it('negative and zero author numbers, insertion codes, a label given twice naming its first position (EC-1 to EC-8)', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'ACDEFGHIK' });
    const labels = ['-3', '-2', '-1', '0', '1', '1A', '1B', '5', '5'];
    b.data.complex('c', {}, (c) => c.unit('A', { of: s }, (u) => u.numbering('auth', { method: 'labels', labels })));
    const u = structureFacts(resolveSpec(b.build()), '/data/c')!.units[0];
    const r = u.residues!;
    expect(Array.from(r.authSeqId)).toEqual([-3, -2, -1, 0, 1, 1, 1, 5, 5]);
    expect(r.insCode).toEqual(['', '', '', '', '', 'A', 'B', '', '']);
    expect(Array.from(r.labeled)).toEqual([1, 1, 1, 1, 1, 1, 1, 1, 1]);
    // a repeated label names its first position (EC-4); every label maps back to a position labeled with it
    expect([u.auth!.position('5'), u.auth!.position('-3'), u.auth!.position('0'), u.auth!.position('1b')]).toEqual([
      7, 0, 3, 6,
    ]);
    for (let p = 0; p < labels.length; p++) expect(labels[u.auth!.position(u.auth!.label(p)!)]).toBe(labels[p]);
    expect([u.auth!.position('-4'), u.auth!.position('1C'), u.auth!.label(9)]).toEqual([-1, -1, undefined]);
  });

  it('residue rows: the author number and insertion code split, observed flags (EC-118)', () => {
    const r = a.residues!;
    expect([r.length, r.authSeqId[10], r.insCode[10], r.labeled[10], r.observed[10]]).toEqual([
      EntLength,
      EntStart + 11,
      'A',
      1,
      1,
    ]);
    expect([r.labeled[0], r.observed[0], r.authSeqId[0], r.insCode[0]]).toEqual([0, 0, 0, '']);
    expect(structuredClone(r)).toEqual(r);
  });
});

/** A spec with the source the nodes came from, and the overlay of those nodes over it. */
function overlayOf(source: string, nodes: readonly { key: string; node: Node }[]): StructureModel {
  const b = createSeqSpec();
  b.data.source(source.split('/').pop()!, { format: 'mmcif', uri: 'x.cif' });
  const model = resolveSpec(b.build());
  const map = new Map<string, NodeEntry>();
  for (const n of nodes) map.set(nodeKey(source, n.key), { source, key: n.key, node: n.node, version: 1 });
  return { model, overlay: buildOverlay(model, { nodes: map, nodeVersion: { value: 1 } }) };
}

const noCode = (c: string | undefined) => (c === undefined || c === '?' || c === '.' ? '' : c);

/** Each described chain's facts against the description's own rows, the description read as written. */
function expectFactsAsDescribed(d: StructureDescription, source: string, model: StructureModel) {
  const facts = structureFacts(model, nodeAddress(source, 'complex'))!;
  const byId = new Map<string, UnitFacts>(facts.units.map((u) => [u.id, u]));
  let chains = 0;
  let labeled = 0;
  for (let i = 0; i < d.chains.label_asym_id.length; i++) {
    const e = d.chains.entity[i];
    const type = d.entities.type[e];
    if (type === 'water') continue;
    chains++;
    const op = d.chains.operator?.[i];
    const id = op ? `${d.chains.label_asym_id[i]}-${op}` : d.chains.label_asym_id[i];
    const u = byId.get(id)!;
    expect(u, id).toBeDefined();
    expect([u.chain, u.authChain, u.operator, u.molecule, u.kind]).toEqual([
      d.chains.label_asym_id[i],
      d.chains.auth_asym_id?.[i],
      op || undefined,
      nodeAddress(source, d.entities.id[e]),
      type === 'polymer' ? 'polymer' : type === 'branched' ? 'branched' : 'other',
    ]);
    // the described rows: a position listed twice is left out by the reading, so it is here
    const rows = new Map<number, number[]>();
    for (let k = 0; k < d.residues.chain.length; k++)
      if (d.residues.chain[k] === i) rows.set(d.residues.position[k], [...(rows.get(d.residues.position[k]) ?? []), k]);
    const r = u.residues!;
    for (let p = 0; p < r.length; p++) {
      const listed = rows.get(p);
      if (listed && listed.length > 1) continue;
      const k = listed?.[0];
      const observed = k === undefined ? 0 : d.residues.observed ? d.residues.observed[k] : 1;
      const auth = k === undefined ? undefined : d.residues.auth_seq_id?.[k];
      const label = auth === undefined ? undefined : `${auth}${noCode(d.residues.ins_code?.[k!])}`;
      if (label !== undefined) labeled++;
      expect(
        [r.observed[p], u.auth?.label(p), r.labeled[p] ? `${r.authSeqId[p]}${r.insCode[p]}` : undefined],
        `${id} ${p}`,
      ).toEqual([observed ? 1 : 0, label, label]);
    }
  }
  expect(facts.units.length).toBe(chains);
  expect(labeled).toBeGreaterThan(chains);
}

describe('a description through structureNodes and structureFacts', () => {
  it('the zinc finger: every chain and residue as described', () => {
    const d = zincFinger();
    expectFactsAsDescribed(d, Source, overlayOf(Source, structureNodes(d, { source: Source }).nodes));
  });

  const viaRuntime = async (id: EntryId) => {
    const source = `/data/${id.toLowerCase()}`;
    const opened = await openEntry(id);
    const { description } = mmcifDescription(
      await frameOf(id),
      { water: false, name: 'description', maxUnits: 10_000 },
      source,
    );
    expectFactsAsDescribed(description, source, overlayOf(source, opened.coverage.nodes!));
  };

  it('1EHZ through the mmCIF runtime: the units’ facts are the file’s rows (EC-118)', () => viaRuntime('1EHZ'));

  describe('the larger entries', { tags: ['medium'] }, () => {
    for (const id of ['1FIN', '4HHB', '3SGJ'] as const)
      it(`${id} through the mmCIF runtime: the units’ facts are the file’s rows (EC-118)`, () => viaRuntime(id));
  });
});
