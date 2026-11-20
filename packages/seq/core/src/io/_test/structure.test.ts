/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Seq io (_spec/README.md "Seq io"): a structure description read into the nodes molecules are read as, equal to the
// same nodes written by hand, valid by the language's rules, keyed by the data's identifiers (a chain removed renames
// nothing; row order changes nothing), and a malformed description said, never thrown.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { validateAndNormalize } from '@molstar/core/repr/schema';
import { hashNode, type Node, nodesEqual } from '@molstar/core/repr/tree';
import { checkDataNode, SeqLanguage, SeqLanguageVersion, SeqValidateOptions } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { type StructureDescription, structureNodes, unitId, unitValues } from '../structure.js';
import { Source, zincFinger, zincFingerNodes } from './fixtures.js';

const read = (d: StructureDescription, water = false) => structureNodes(d, { source: Source, water });
const byKey = (d: StructureDescription) => new Map(read(d).nodes.map((n) => [n.key, n.node]));

describe('structureNodes', () => {
  it('reads a description into the nodes written by hand: a sequence, a molecule, compounds, a complex', () => {
    const { nodes, issues } = read(zincFinger());
    const hand = zincFingerNodes();
    expect(nodes.map((n) => n.key)).toEqual(['1', '2', 'compounds', 'complex']);
    expect(nodes.map((n) => n.node.kind)).toEqual(hand.map((n) => n.kind));
    for (let i = 0; i < hand.length; i++) expect(nodes[i].node).toEqual(hand[i]);
    expect(issues.map((i) => [i.code, i.message])).toEqual([
      ['seq.io.structure.letters', '1 positions’ components have no standard letter: X or N there'],
    ]);
  });

  it('gives nodes the language takes under data and its rules pass, as the overlay checks a host’s', () => {
    for (const { node } of read(zincFinger(), true).nodes) {
      const spec: Node = {
        kind: 'spec',
        params: { version: SeqLanguageVersion },
        children: [{ kind: 'data', id: 'data', children: [node] }],
      };
      expect(validateAndNormalize(SeqLanguage, spec, SeqValidateOptions).issues).toEqual([]);
      const said: string[] = [];
      checkDataNode(node, { target: () => undefined }, (code) => said.push(code));
      expect(said).toEqual([]);
    }
  });

  it('reads water only when asked', () => {
    expect(read(zincFinger()).nodes.map((n) => n.key)).not.toContain('3');
    const withWater = read(zincFinger(), true);
    expect(withWater.nodes.map((n) => n.key)).toContain('3');
    expect(withWater.systems.map((s) => s.name)).toEqual(['A-1_555', 'B-1_555', 'C-1_555', 'D-1_555']);
  });

  it('keeps every other key and node when a chain is removed (a host’s edit renames nothing, EC-181)', () => {
    const d = zincFinger();
    const before = byKey(d);
    const keep = (i: number) => d.residues.chain[i] !== 1;
    const rows = [...d.residues.chain.keys()].filter(keep);
    const remap = (c: number) => (c > 1 ? c - 1 : c);
    const removed: StructureDescription = {
      ...d,
      chains: {
        label_asym_id: d.chains.label_asym_id.filter((_, i) => i !== 1),
        auth_asym_id: d.chains.auth_asym_id!.filter((_, i) => i !== 1),
        entity: d.chains.entity.filter((_, i) => i !== 1),
        operator: d.chains.operator!.filter((_, i) => i !== 1),
      },
      residues: {
        chain: Int32Array.from(rows, (i) => remap(d.residues.chain[i])),
        position: Int32Array.from(rows, (i) => d.residues.position[i]),
        auth_seq_id: Int32Array.from(rows, (i) => d.residues.auth_seq_id![i]),
        ins_code: rows.map((i) => d.residues.ins_code![i]),
      },
      links: { ...d.links!, a_chain: Int32Array.of(0), b_chain: Int32Array.of(1) },
      values: undefined,
    };
    const after = byKey(removed);
    expect([...after.keys()]).toEqual([...before.keys()]);
    for (const key of ['1', '2', 'compounds']) expect(nodesEqual(after.get(key)!, before.get(key)!)).toBe(true);
    const units = (n: Node) => n.children!.filter((c) => c.kind === 'unit');
    expect(units(after.get('complex')!).map((u) => u.id)).toEqual(['A-1_555', 'C-1_555']);
    for (const u of units(after.get('complex')!))
      expect(hashNode(u)).toBe(hashNode(units(before.get('complex')!).find((b) => b.id === u.id)!));
  });

  it('a property: row order changes no key and no node (seeded permutations of every table)', () => {
    const d = zincFinger();
    const base = byKey(d);
    const rng = new PCG(7);
    const permutation = (n: number) => {
      const p = Array.from({ length: n }, (_, i) => i);
      for (let i = n - 1; i > 0; i--) {
        const j = rng.int() % (i + 1);
        [p[i], p[j]] = [p[j], p[i]];
      }
      return p;
    };
    for (let run = 0; run < 25; run++) {
      const e = permutation(d.entities.id.length);
      const inv = e.map((_, i) => e.indexOf(i));
      const c = permutation(d.components.entity.length);
      const r = permutation(d.residues.chain.length);
      const k = permutation(d.compounds!.id.length);
      const shuffled: StructureDescription = {
        ...d,
        entities: {
          id: e.map((i) => d.entities.id[i]),
          type: e.map((i) => d.entities.type[i]),
          polymerType: e.map((i) => d.entities.polymerType![i]),
          name: e.map((i) => d.entities.name![i]),
          letters: e.map((i) => d.entities.letters![i]),
        },
        components: {
          entity: Int32Array.from(c, (i) => inv[d.components.entity[i]]),
          position: Int32Array.from(c, (i) => d.components.position[i]),
          component: c.map((i) => d.components.component[i]),
        },
        identifiers: { ...d.identifiers!, entity: d.identifiers!.entity.map((x) => inv[x]) },
        chains: { ...d.chains, entity: d.chains.entity.map((x) => inv[x]) },
        residues: {
          chain: Int32Array.from(r, (i) => d.residues.chain[i]),
          position: Int32Array.from(r, (i) => d.residues.position[i]),
          auth_seq_id: Int32Array.from(r, (i) => d.residues.auth_seq_id![i]),
          ins_code: r.map((i) => d.residues.ins_code![i]),
        },
        compounds: {
          id: k.map((i) => d.compounds!.id[i]),
          name: k.map((i) => d.compounds!.name![i]),
          parents: k.map((i) => d.compounds!.parents![i]),
        },
        atoms: { ...d.atoms!, entity: d.atoms!.entity.map((x) => inv[x]) },
      };
      const got = byKey(shuffled);
      expect([...got.keys()]).toEqual([...base.keys()]);
      for (const [key, node] of got) expect(hashNode(node), `run ${run}, ${key}`).toBe(hashNode(base.get(key)!));
    }
  });

  it('says a malformed description and leaves out what it cannot read, never throwing', () => {
    const d = zincFinger();
    const bad: StructureDescription = {
      ...d,
      entities: { ...d.entities, id: ['1', 'two words', '3'] },
      components: { ...d.components, entity: Int32Array.of(0, 0, 0, 0, 0, 9, 2) },
      chains: { ...d.chains, entity: Int32Array.of(0, 0, 7, 2) },
      residues: { ...d.residues, ins_code: ['', ''] },
      links: { ...d.links!, a_position: Int32Array.of(99) },
    };
    const { nodes, issues } = read(bad);
    expect(nodes.map((n) => n.key)).toEqual(['1', 'compounds', 'complex']);
    expect(issues.map((i) => i.code).sort()).toEqual([
      'seq.io.structure.columns',
      'seq.io.structure.id',
      'seq.io.structure.letters',
      'seq.io.structure.rows',
      'seq.io.structure.rows',
    ]);
    expect(issues.every((i) => i.address === Source)).toBe(true);
  });

  it('leaves a molecule’s hydrogens out, said (EC-64)', () => {
    const d = zincFinger();
    const withH: StructureDescription = {
      ...d,
      atoms: { entity: Int32Array.of(1, 1), residue: Int32Array.of(0, 0), name: ['ZN', 'H1'], element: ['ZN', 'H'] },
    };
    const { nodes, issues } = read(withH);
    expect(hashNode(nodes[1].node)).toBe(hashNode(zincFingerNodes()[1]));
    expect(issues.map((i) => i.code)).toContain('seq.io.structure.hydrogen');
  });
});

describe('unitId and unitValues', () => {
  it('a unit’s id is its chain and operator, what an id takes', () => {
    expect(unitId('A', undefined)).toBe('A');
    expect(unitId('A', '1_555')).toBe('A-1_555');
    expect(unitId('A B', 'X/1')).toBe('A_B-X_1');
  });

  it('a named value over each unit’s positions, NaN where a residue has none', () => {
    const d = zincFinger();
    const v = unitValues(d, 'b_factor', read(d).systems)!;
    expect([...v.keys()]).toEqual(['A-1_555', 'B-1_555', 'C-1_555']);
    expect(Array.from(v.get('A-1_555')!)).toEqual([Number.NaN, 20, 21, 22, 23]);
    expect(unitValues(d, 'plddt', read(d).systems)).toBeUndefined();
  });
});
