/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A source's molecules in the overlay: the zinc finger's structure description read as nodes (`structureNodes`), its
// sequence, ion and complex held as systems, units checked by the language and bounded by `maxUnits`, reached by the
// query forms as a spec's are.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { type Node, printAddress } from '@molstar/core/repr/tree';
import { standardLetter } from '@molstar/core/util/sequence-constants';
import { type NodeEntry, nodeKey } from '@molstar/seqx-core/data';
import { structureNodes } from '@molstar/seqx-core/io/structure';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { evaluateQuery, partOf } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { Source, zincFinger } from '../../io/_test/fixtures.js';
import { resolveSpec } from '../model.js';
import { modifiedIndex } from '../molecules.js';
import { buildOverlay, nodeAddress } from '../overlay.js';
import { overlayQueryModel } from '../query-model.js';

const Complex = nodeAddress(Source, 'complex');
const UnitA = `${Complex}/A-1_555`;
const UnitB = `${Complex}/B-1_555`;
const Seq1 = nodeAddress(Source, '1');

function spec(): Node {
  const b = createSeqSpec();
  b.data.source('zf', { format: 'mmcif', uri: 'zf.cif' });
  return b.build();
}

function reader(nodes: readonly { key: string; node: Node }[]) {
  const map = new Map<string, NodeEntry>();
  for (const n of nodes) map.set(nodeKey(Source, n.key), { source: Source, key: n.key, node: n.node, version: 1 });
  return { nodes: map, nodeVersion: { value: 1 } };
}

function setup(nodes = structureNodes(zincFinger(), { source: Source }).nodes, maxUnits?: number) {
  const s = spec();
  const model = resolveSpec(s, maxUnits ? { config: { maxUnits } } : {});
  const overlay = buildOverlay(model, reader(nodes));
  const query = overlayQueryModel(model, overlay);
  const run = (address: string) => evaluateQuery(s, query, address, { limit: 100 });
  return { model, overlay, query, run };
}

const loc = (r: ReturnType<ReturnType<typeof setup>['run']>, system: string) =>
  r.result?.type === 'location' ? Array.from(partOf(r.result.location, system)?.ranges ?? []) : r.issues;

const nodes = (r: ReturnType<ReturnType<typeof setup>['run']>) =>
  r.result?.type === 'nodes' ? r.result.nodes.map((m) => printAddress(m.address)) : r.issues;

describe('the overlay’s molecules', () => {
  it('a source’s sequence, molecule and complex are systems: units on their molecule from 0 (EC-35)', () => {
    const { overlay } = setup();
    expect(overlay.issues.filter((i) => i.severity === 'error')).toEqual([]);
    expect(overlay.system(Seq1)).toMatchObject({ kind: 'sequence', length: 5, letters: 'CHMCX' });
    expect(overlay.system(nodeAddress(Source, '2'))).toMatchObject({ kind: 'molecule', length: 1 });
    expect(overlay.system(`${nodeAddress(Source, '2')}/(atoms)`)).toMatchObject({ kind: 'atoms', length: 1 });
    expect(overlay.system(UnitA)).toMatchObject({ kind: 'unit', length: 5, parent: { system: Seq1, start: 0 } });
    const c = overlay.molecular.complexes.get(Complex)!;
    expect([...c.copies.values()]).toEqual([[UnitA, UnitB], [`${Complex}/C-1_555`]]);
    const links = overlay.annotation(`${Complex}/links`)!;
    expect([links.units?.length, Array.from(links.pairs!.aUnit!), Array.from(links.pairs!.bUnit!)]).toEqual([
      3,
      [0, 1],
      [2, 2],
    ]);
  });

  it('paths: a unit steps to its molecule; no path joins two copies (EC-25)', () => {
    const { overlay } = setup();
    expect(overlay.path(UnitA, Seq1).path?.edges.map((e) => e.unit)).toEqual(['up']);
    expect(overlay.path(UnitA, UnitB).issue?.code).toBe('seq.resolve.path.none');
  });

  it('what a unit is drawn through is held with it: its complex and the molecules its units instantiate', () => {
    const { overlay } = setup();
    expect(
      overlay
        .holds([UnitA])
        .map((n) => n.key)
        .sort(),
    ).toEqual(['1', '2', 'complex']);
    expect(overlay.holds([`${nodeAddress(Source, '2')}/(atoms)`]).map((n) => n.key)).toEqual(['2']);
  });

  it('a unit’s auth numbering resolves through the chain form: (chain :id A)/(at :label 12A :numbering auth)', () => {
    const { run } = setup();
    expect(nodes(run(`${Complex}/(chain :id A)`))).toEqual([UnitA]);
    expect(loc(run(`${Complex}/(chain :id A)/(at :label 12A :numbering auth)`), UnitA)).toEqual([3, 4]);
    expect(nodes(run(`${Complex}/(chain :auth A :operator 1_555)`))).toEqual([UnitA, `${Complex}/C-1_555`]);
    expect(run(`${Complex}/(chain :id Q)`).issues.map((i) => i.code)).toEqual(['seq.query.node.missing']);
  });

  it('components, entities and molecules of a source’s complex', () => {
    const { run } = setup();
    expect(loc(run(`${UnitA}/(elements :component MSE)`), UnitA)).toEqual([2, 3]);
    expect(loc(run(`${Seq1}/(elements :standard false)`), Seq1)).toEqual([2, 3, 4, 5]);
    expect(nodes(run(`${Complex}/(molecule :entity 2)`))).toEqual([nodeAddress(Source, '2')]);
    expect(nodes(run(`${Complex}/(molecule :identifier "ccd:ZN")`))).toEqual([nodeAddress(Source, '2')]);
    expect(loc(run(`${Complex}/C-1_555/(elements :atom ZN)`), `${Complex}/C-1_555/(atoms)`)).toEqual([0, 1]);
  });

  it('a complex past maxUnits is refused whole, with an issue', () => {
    const { overlay } = setup(undefined, 2);
    expect(overlay.system(UnitA)).toBeUndefined();
    expect(overlay.issues.map((i) => [i.code, i.address])).toEqual([['seq.resolve.complex.units', Complex]]);
  });

  it('a unit or links failing the language’s checks are left out with their issues (checkUnit, checkLinks)', () => {
    const complex: Node = {
      kind: 'complex',
      id: 'complex',
      children: [
        { kind: 'unit', id: 'A', params: { of: Seq1, range: [4, 2] } },
        { kind: 'unit', id: 'B', params: { of: Seq1, range: [0, 5] } },
        {
          kind: 'pairs',
          id: 'links',
          params: {
            rows: { a_unit: ['B'], a_start: [0], a_end: [1], b_unit: ['Z'], b_start: [0], b_end: [1] },
          },
        },
      ],
    };
    const seq = structureNodes(zincFinger(), { source: Source }).nodes.find((n) => n.key === '1')!;
    const { overlay } = setup([seq, { key: 'complex', node: complex }]);
    expect(overlay.system(`${Complex}/A`)).toBeUndefined();
    expect(overlay.system(`${Complex}/B`)).toMatchObject({ length: 5 });
    expect(overlay.annotation(`${Complex}/links`)).toBeUndefined();
    expect(overlay.issues.map((i) => [i.code, i.address])).toEqual([
      ['seq.lang.unit.range', `${Complex}/A`],
      ['seq.lang.link.unit', `${Complex}/links`],
    ]);
  });
});

describe('properties of a source’s modified positions', () => {
  // components a position: standard ones, and modified ones whose compound names one parent
  const Modified: Record<string, string> = { MSE: 'MET', SEP: 'SER', TPO: 'THR', PTR: 'TYR', CSO: 'CYS' };
  const Standard = ['ALA', 'GLY', 'MET', 'SER', 'THR', 'TYR', 'CYS', 'LYS'];

  it('a modified position’s letter is its parent’s', () => {
    const rng = new PCG(104);
    for (let k = 0; k < 40; k++) {
      const n = 1 + (rng.int() % 12);
      const component = Array.from({ length: n }, () =>
        rng.int() % 3 ? Standard[rng.int() % Standard.length] : Object.keys(Modified)[rng.int() % 5],
      );
      const d = zincFinger();
      const read = structureNodes(
        {
          ...d,
          entities: { id: ['1'], type: ['polymer'], polymerType: ['polypeptide(L)'] },
          components: {
            entity: new Int32Array(n),
            position: Int32Array.from({ length: n }, (_, i) => i),
            component,
          },
          identifiers: undefined,
          chains: { label_asym_id: ['A'], entity: Int32Array.of(0) },
          residues: { chain: new Int32Array(n), position: Int32Array.from({ length: n }, (_, i) => i) },
          links: undefined,
          compounds: { id: Object.keys(Modified), parents: Object.values(Modified).map((p) => [p]) },
          atoms: undefined,
          values: undefined,
        },
        { source: Source },
      );
      const node = read.nodes.find((x) => x.key === '1')!.node;
      const letters = node.params!.letters as string;
      const index = modifiedIndex(node);
      for (let i = 0; i < (index?.position.length ?? 0); i++) {
        const parent = Modified[index!.component[i]];
        expect(letters[index!.position[i]]).toBe(standardLetter('protein', parent));
      }
      expect(index?.position.length ?? 0).toBe(component.filter((c) => c in Modified).length);
    }
  });
});
