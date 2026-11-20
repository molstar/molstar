/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, identifiers, modified } from '@molstar/seqx-core/lang';
import type { StructureDescription } from '../structure.js';

// A small structure as a host describes it, and the same nodes written by hand with the language's builder: a zinc
// finger's two chains (MSE at 2, residue 0 unmodeled, an insertion code), a zinc, a water, two metal links (rows out of order, as
// identifiers).

export const Source = '/data/zf';

export function zincFinger(): StructureDescription {
  return {
    name: 'Zinc finger',
    entities: {
      id: ['1', '2', '3'],
      type: ['polymer', 'non-polymer', 'water'],
      polymerType: ['polypeptide(L)', undefined, undefined],
      name: ['Finger', 'ZINC ION', 'water'],
      letters: [undefined, undefined, undefined],
    },
    components: {
      entity: Int32Array.of(0, 0, 0, 0, 0, 1, 2),
      position: Int32Array.of(0, 1, 2, 3, 4, 0, 0),
      component: ['CYS', 'HIS', 'MSE', 'CYS', 'XYZ', 'ZN', 'HOH'],
    },
    identifiers: { entity: Int32Array.of(1, 1), namespace: ['pubchem', 'ccd'], value: ['32051', 'ZN'] },
    chains: {
      label_asym_id: ['A', 'B', 'C', 'D'],
      auth_asym_id: ['A', 'B', 'A', 'A'],
      entity: Int32Array.of(0, 0, 1, 2),
      operator: ['1_555', '1_555', '1_555', '1_555'],
    },
    residues: {
      chain: Int32Array.of(0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 3),
      position: Int32Array.of(1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 0),
      auth_seq_id: Int32Array.of(11, 12, 12, 13, 10, 11, 12, 13, 14, 101, 201),
      ins_code: ['', '', 'A', '', '', '', '', '', '', '', ''],
    },
    links: {
      a_chain: Int32Array.of(1, 0),
      a_position: Int32Array.of(3, 3),
      b_chain: Int32Array.of(2, 2),
      b_position: Int32Array.of(0, 0),
      type: ['metal', 'metal'],
      a_atom: ['SG', 'SG'],
      b_atom: ['ZN', 'ZN'],
    },
    compounds: { id: ['ZN', 'MSE'], name: ['ZINC ION', 'SELENOMETHIONINE'], parents: [undefined, ['MET']] },
    atoms: { entity: Int32Array.of(1), residue: Int32Array.of(0), name: ['ZN'], element: ['ZN'] },
    values: { b_factor: Float32Array.of(20, 21, 22, 23, 30, 31, 32, 33, 34, 50, 60) },
  };
}

/** The nodes `zincFinger` is read as, written with the builder: what the default adapter and every host give. */
export function zincFingerNodes(): readonly Node[] {
  const b = createSeqSpec();
  b.data.sequence('1', {
    alphabet: 'protein',
    letters: 'CHMCX',
    modified: modified([
      { position: 2, component: 'MSE' },
      { position: 4, component: 'XYZ' },
    ]),
  });
  b.data.molecule('2', {
    class: 'ion',
    name: 'ZINC ION',
    residues: ['ZN'],
    atoms: { residue: [0], name: ['ZN'], element: ['ZN'] },
    identifiers: identifiers({ ccd: 'ZN', pubchem: '32051' }),
  });
  b.data.compound('MSE', { name: 'SELENOMETHIONINE', parents: ['MET'] });
  b.data.compound('ZN', { name: 'ZINC ION' });
  b.data.complex('complex', { name: 'Zinc finger' }, (c) => {
    const auth = (labels: (string | null)[]) => ({ method: 'labels' as const, scheme: 'auth', labels });
    c.unit(
      'A-1_555',
      {
        of: `${Source}/(node :key 1)`,
        chain: 'A',
        auth_chain: 'A',
        operator: '1_555',
        observed: { start: [1], end: [5] },
      },
      (u) => u.numbering('auth', auth([null, '11', '12', '12A', '13'])),
    );
    c.unit('B-1_555', { of: `${Source}/(node :key 1)`, chain: 'B', auth_chain: 'B', operator: '1_555' }, (u) =>
      u.numbering('auth', auth(['10', '11', '12', '13', '14'])),
    );
    c.unit('C-1_555', { of: `${Source}/(node :key 2)`, chain: 'C', auth_chain: 'A', operator: '1_555' }, (u) =>
      u.numbering('auth', auth(['101'])),
    );
    c.links('links', {
      rows: {
        a_unit: ['A-1_555', 'B-1_555'],
        a_start: [3, 3],
        a_end: [4, 4],
        b_unit: ['C-1_555', 'C-1_555'],
        b_start: [0, 0],
        b_end: [1, 1],
        type: ['metal', 'metal'],
        a_atom: ['SG', 'SG'],
        b_atom: ['ZN', 'ZN'],
      },
    });
  });
  return b.build().children?.[0].children ?? [];
}

/** A description without one chain row: its residues, values and links with it, the other rows' indices moved. */
export function withoutChain(d: StructureDescription, chain: number): StructureDescription {
  const keep = <T>(xs: ArrayLike<T>, ok: (i: number) => boolean) => Array.from(xs).filter((_, i) => ok(i));
  const moved = (c: number) => (c > chain ? c - 1 : c);
  const rows = (i: number) => d.residues.chain[i] !== chain;
  const l = d.links!;
  const links = (i: number) => l.a_chain[i] !== chain && l.b_chain[i] !== chain;
  const c = d.chains;
  return {
    ...d,
    chains: {
      label_asym_id: keep(c.label_asym_id, (i) => i !== chain),
      auth_asym_id: keep(c.auth_asym_id!, (i) => i !== chain),
      entity: Int32Array.from(keep(c.entity, (i) => i !== chain)),
      operator: keep(c.operator!, (i) => i !== chain),
    },
    residues: {
      chain: Int32Array.from(keep(d.residues.chain, rows), moved),
      position: Int32Array.from(keep(d.residues.position, rows)),
      auth_seq_id: Int32Array.from(keep(d.residues.auth_seq_id!, rows)),
      ins_code: keep(d.residues.ins_code!, rows),
    },
    links: {
      a_chain: Int32Array.from(keep(l.a_chain, links), moved),
      a_position: Int32Array.from(keep(l.a_position, links)),
      b_chain: Int32Array.from(keep(l.b_chain, links), moved),
      b_position: Int32Array.from(keep(l.b_position, links)),
      type: keep(l.type!, links),
      a_atom: keep(l.a_atom!, links),
      b_atom: keep(l.b_atom!, links),
    },
    values: Object.fromEntries(Object.entries(d.values ?? {}).map(([k, v]) => [k, Float32Array.from(keep(v, rows))])),
  };
}

/**
 * The same description with every table's rows in another order (a seeded permutation each), the indices between them
 * moved with them; a position's alternatives keep their order (it is data: the first gives the letter).
 */
export function permuted(d: StructureDescription, next: () => number): StructureDescription {
  const perm = (n: number) => {
    const p = Array.from({ length: n }, (_, i) => i);
    for (let i = n - 1; i > 0; i--) {
      const j = next() % (i + 1);
      [p[i], p[j]] = [p[j], p[i]];
    }
    return p;
  };
  const inverse = (p: number[]) => {
    const out = new Array<number>(p.length);
    p.forEach((from, to) => (out[from] = to));
    return out;
  };
  const pick = <T>(xs: ArrayLike<T> | undefined, p: number[]) => (xs ? p.map((i) => xs[i]) : undefined);
  const ints = (xs: ArrayLike<number>, p: number[], map = (x: number) => x) => Int32Array.from(p, (i) => map(xs[i]));
  const e = perm(d.entities.id.length);
  const eTo = inverse(e);
  const ch = perm(d.chains.label_asym_id.length);
  const chTo = inverse(ch);
  const c = perm(d.components.entity.length);
  const r = perm(d.residues.chain.length);
  const id = perm(d.identifiers!.entity.length);
  const l = perm(d.links!.a_chain.length);
  const k = perm(d.compounds!.id.length);
  const a = perm(d.atoms!.entity.length);
  const E = (x: number) => eTo[x];
  const C = (x: number) => chTo[x];
  return {
    ...d,
    entities: {
      id: pick(d.entities.id, e)!,
      type: pick(d.entities.type, e)!,
      polymerType: pick(d.entities.polymerType, e),
      name: pick(d.entities.name, e),
      letters: pick(d.entities.letters, e),
    },
    components: {
      entity: ints(d.components.entity, c, E),
      position: ints(d.components.position, c),
      component: pick(d.components.component, c)!,
    },
    identifiers: {
      entity: ints(d.identifiers!.entity, id, E),
      namespace: pick(d.identifiers!.namespace, id)!,
      value: pick(d.identifiers!.value, id)!,
    },
    chains: {
      label_asym_id: pick(d.chains.label_asym_id, ch)!,
      auth_asym_id: pick(d.chains.auth_asym_id, ch),
      entity: ints(d.chains.entity, ch, E),
      operator: pick(d.chains.operator, ch),
    },
    residues: {
      chain: ints(d.residues.chain, r, C),
      position: ints(d.residues.position, r),
      auth_seq_id: ints(d.residues.auth_seq_id!, r),
      ins_code: pick(d.residues.ins_code, r),
    },
    links: {
      a_chain: ints(d.links!.a_chain, l, C),
      a_position: ints(d.links!.a_position, l),
      b_chain: ints(d.links!.b_chain, l, C),
      b_position: ints(d.links!.b_position, l),
      type: pick(d.links!.type, l),
      a_atom: pick(d.links!.a_atom, l),
      b_atom: pick(d.links!.b_atom, l),
    },
    compounds: {
      id: pick(d.compounds!.id, k)!,
      name: pick(d.compounds!.name, k),
      parents: pick(d.compounds!.parents, k),
    },
    atoms: {
      entity: ints(d.atoms!.entity, a, E),
      residue: ints(d.atoms!.residue, a),
      name: pick(d.atoms!.name, a)!,
      element: pick(d.atoms!.element, a)!,
    },
    values: Object.fromEntries(Object.entries(d.values ?? {}).map(([n, v]) => [n, Float32Array.from(r, (i) => v[i])])),
  };
}
