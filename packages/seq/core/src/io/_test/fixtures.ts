/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, identifiers, modified } from '@molstar/seqx-core/lang';
import type { StructureDescription } from '../structure.js';

// A small structure as a host describes it, and the same nodes written by hand with the language's builder: a zinc
// finger's two chains (MSE at 2, residue 0 unmodeled, an insertion code), a zinc, waters, a metal link.

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
    identifiers: { entity: Int32Array.of(1), namespace: ['ccd'], value: ['ZN'] },
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
      a_chain: Int32Array.of(0),
      a_position: Int32Array.of(3),
      b_chain: Int32Array.of(2),
      b_position: Int32Array.of(0),
      type: ['metal'],
      a_atom: ['SG'],
      b_atom: ['ZN'],
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
    identifiers: identifiers({ ccd: 'ZN' }),
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
        a_unit: ['A-1_555'],
        a_start: [3],
        a_end: [4],
        b_unit: ['C-1_555'],
        b_start: [0],
        b_end: [1],
        type: ['metal'],
        a_atom: ['SG'],
        b_atom: ['ZN'],
      },
    });
  });
  return b.build().children?.[0].children ?? [];
}
