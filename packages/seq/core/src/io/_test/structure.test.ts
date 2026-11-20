/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Seq io (_spec/README.md "Seq io"): a structure description read into the nodes molecules are read as, equal to the
// same nodes written by hand, valid by the language's rules, keyed and ordered by the data's identifiers (a chain
// removed renames nothing; row order changes nothing), and a malformed or oversized description said within its
// bounds, never thrown.

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { validateAndNormalize } from '@molstar/core/repr/schema';
import { hashNode, type Node, nodesEqual } from '@molstar/core/repr/tree';
import { checkDataNode, SeqLanguage, SeqLanguageVersion, SeqValidateOptions } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { DefaultStructureLimits, type StructureDescription, structureNodes } from '../structure.js';
import { permuted, Source, withoutChain, zincFinger, zincFingerNodes } from './fixtures.js';

const read = (d: StructureDescription, water = false) => structureNodes(d, { source: Source, water });
const byKey = (d: StructureDescription) => new Map(read(d).nodes.map((n) => [n.key, n.node]));
const units = (n: Node | undefined) => (n?.children ?? []).filter((c) => c.kind === 'unit');

/** The nodes pass the language's description and its data rules, as the overlay checks a host's. */
function expectValid(nodes: readonly { node: Node }[]) {
  for (const { node } of nodes) {
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
}

describe('structureNodes', () => {
  it('reads a description into the nodes written by hand: a sequence, a molecule, compounds, a complex', () => {
    const { nodes, issues } = read(zincFinger());
    const hand = zincFingerNodes();
    expect(nodes.map((n) => n.key)).toEqual(['1', '2', 'compounds', 'complex']);
    for (let i = 0; i < hand.length; i++) expect(nodes[i].node).toEqual(hand[i]);
    expect(issues.map((i) => [i.code, i.message])).toEqual([
      ['seq.io.structure.letters', '1 positions’ components have no standard letter: X or N there'],
    ]);
  });

  it('gives nodes the language takes under data and its rules pass, as the overlay checks a host’s', () => {
    expectValid(read(zincFinger(), true).nodes);
  });

  it('reads water only when asked', () => {
    expect(read(zincFinger()).nodes.map((n) => n.key)).not.toContain('3');
    const withWater = read(zincFinger(), true);
    expect(withWater.nodes.map((n) => n.key)).toContain('3');
    expect(withWater.systems.map((s) => s.name)).toEqual(['A-1_555', 'B-1_555', 'C-1_555', 'D-1_555']);
  });

  it('keeps every other key and node when a chain is removed (a host’s edit renames nothing, EC-181)', () => {
    const before = byKey(zincFinger());
    const after = byKey(withoutChain(zincFinger(), 1));
    expect([...after.keys()]).toEqual([...before.keys()]);
    for (const key of ['1', '2', 'compounds']) expect(nodesEqual(after.get(key)!, before.get(key)!)).toBe(true);
    expect(units(after.get('complex')).map((u) => u.id)).toEqual(['A-1_555', 'C-1_555']);
    for (const u of units(after.get('complex')))
      expect(hashNode(u)).toBe(hashNode(units(before.get('complex')).find((b) => b.id === u.id)!));
  });

  it('a property: row order changes no key and no content (seeded permutations of every table, chains, links, atoms and identifiers too)', () => {
    const d = zincFinger();
    const base = byKey(d);
    const baseValues = read(d).values('b_factor');
    const rng = new PCG(7);
    for (let run = 0; run < 25; run++) {
      const shuffled = permuted(d, () => rng.int());
      const got = byKey(shuffled);
      expect([...got.keys()]).toEqual([...base.keys()]);
      for (const [key, node] of got) expect(hashNode(node), `run ${run}, ${key}`).toBe(hashNode(base.get(key)!));
      expect(read(shuffled).values('b_factor')).toEqual(baseValues);
    }
  });

  it('says a malformed description and leaves out what it cannot read, never throwing', () => {
    const d = zincFinger();
    const bad: StructureDescription = {
      ...d,
      entities: { ...d.entities, id: ['1', 'two words', '3'] },
      components: { ...d.components, entity: Int32Array.of(0, 0, 0, 0, 0, 9, 2) },
      residues: { ...d.residues, ins_code: ['', ''] },
      links: { ...d.links!, a_position: Int32Array.of(3, 99) },
    };
    const { nodes, issues } = read(bad);
    expect(nodes.map((n) => n.key)).toEqual(['1', 'compounds', 'complex']);
    expect(issues.map((i) => i.code).sort()).toEqual([
      'seq.io.structure.columns',
      'seq.io.structure.id',
      'seq.io.structure.letters',
      'seq.io.structure.rows',
    ]);
    expect(issues.every((i) => i.address === Source)).toBe(true);
  });

  it('a molecule’s own links are its links child, by their ends; a polymer’s or a position past the molecule left out, said', () => {
    const d = zincFinger();
    const withLinks: StructureDescription = {
      ...d,
      entityLinks: {
        entity: Int32Array.of(1, 1, 0, 1),
        a_position: Int32Array.of(0, 0, 0, 0),
        b_position: Int32Array.of(0, 0, 1, 5),
        type: ['glycosidic', 'covalent', 'covalent', 'covalent'],
        a_atom: ['O4', null, null, null],
        b_atom: ['C1', null, null, null],
        order: ['sing', null, null, null],
      },
    };
    const { nodes, issues } = read(withLinks);
    const links = nodes[1].node.children?.find((c) => c.kind === 'pairs');
    const rows = links?.params?.rows as Record<string, any>;
    expect(links?.id).toBe('links');
    expect(Array.from(rows.a_start)).toEqual([0, 0]);
    expect([rows.type.get(0), rows.type.get(1), rows.a_atom.get(1), rows.order.get(1)]).toEqual([
      'covalent',
      'glycosidic',
      'O4',
      'sing',
    ]);
    expect(nodes[0].node.children).toBeUndefined();
    expect(issues.map((i) => i.code)).toContain('seq.io.structure.rows');
    expectValid(nodes);
  });

  it('a non-polymer is an ion by its component’s formula, not by the atoms modeled (a phosphate with its P alone)', () => {
    const d = zincFinger();
    const phosphate: StructureDescription = {
      ...d,
      components: { ...d.components, component: ['CYS', 'HIS', 'MSE', 'CYS', 'XYZ', 'PO4', 'HOH'] },
      compounds: { id: ['PO4', 'ZN'], formula: ['O4 P -3', 'Zn 2'] },
      atoms: { entity: Int32Array.of(1), residue: Int32Array.of(0), name: ['P'], element: ['P'] },
    };
    expect(read(phosphate).nodes[1].node.params?.class).toBe('ligand');
    const zinc: StructureDescription = { ...d, compounds: { id: ['ZN'], formula: ['Zn 2'] }, atoms: undefined };
    expect(read(zinc).nodes[1].node.params?.class).toBe('ion');
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

describe('bounds', () => {
  it('a position past the bound is said and left out, never sized by (no array of two billion)', () => {
    const d = zincFinger();
    const far: StructureDescription = {
      ...d,
      components: { ...d.components, position: Int32Array.of(0, 1, 2, 3, 2e9, 0, 0) },
    };
    const { nodes, issues } = read(far);
    expect((nodes[0].node.params?.letters as string).length).toBe(4);
    expect(issues.map((i) => i.code)).toContain('seq.io.structure.bounds');
  });

  it('units, rows and lengths past the configured limits are said and left out', () => {
    const d = zincFinger();
    const few = structureNodes(d, { source: Source, limits: { maxUnits: 2 } });
    expect(few.systems.map((s) => s.name)).toEqual(['A-1_555', 'B-1_555']);
    expect(few.issues.map((i) => i.code)).toContain('seq.io.structure.bounds');
    // components (7 rows) and residues (11) past 5 rows: left out, the polymer read with no positions
    const rows = structureNodes(d, { source: Source, limits: { maxRows: 5 } });
    expect(rows.nodes[0].node.params?.letters).toBe('');
    expect(rows.issues.filter((i) => i.code === 'seq.io.structure.bounds').length).toBe(2);
    const short = structureNodes(d, { source: Source, limits: { maxLength: 3 } });
    expect((short.nodes[0].node.params?.letters as string).length).toBe(3);
  });

  it('a description with tables missing or of other types reads as empty, said', () => {
    for (const broken of [
      {},
      { entities: 'none', chains: null },
      { ...zincFinger(), residues: undefined },
      { ...zincFinger(), chains: { label_asym_id: 'A', entity: [0] } },
      null,
    ]) {
      const got = read(broken as unknown as StructureDescription);
      expect(got.issues.map((i) => i.code)).toContain('seq.io.structure.columns');
      expectValid(got.nodes);
    }
  });

  it('a fuzz: seeded malformed descriptions never throw, say only coded issues, and give valid nodes', () => {
    const rng = new PCG(11);
    const junk = [undefined, null, -1, 2e9, 1.5, Number.NaN, '', 'two words', 'A-1', {}, [], [1, 2], 'H', 0, 7];
    const pick = <T>(xs: readonly T[]) => xs[rng.int() % xs.length];
    for (let run = 0; run < 200; run++) {
      const d = JSON.parse(
        JSON.stringify(zincFinger(), (_, v) => (ArrayBuffer.isView(v) ? Array.from(v as never) : v)),
      );
      for (let n = 1 + (rng.int() % 4); n > 0; n--) {
        const table = pick(Object.keys(d));
        const t = d[table];
        if (typeof t !== 'object' || t === null) d[table] = pick(junk);
        else {
          const column = pick(Object.keys(t));
          if (rng.int() % 3 === 0 || !Array.isArray(t[column])) t[column] = pick(junk);
          else t[column][rng.int() % Math.max(1, t[column].length)] = pick(junk);
        }
      }
      let got: ReturnType<typeof read> | undefined;
      expect(() => (got = read(d, rng.int() % 2 === 0)), `run ${run}`).not.toThrow();
      expect(got!.issues.every((i) => i.code.startsWith('seq.io.structure.'))).toBe(true);
      expectValid(got!.nodes);
      for (const s of got!.systems) expect(s.length).toBeLessThanOrEqual(DefaultStructureLimits.maxLength);
    }
  });
});

describe('units and their values', () => {
  it('a unit’s id is its chain and its operator; one that makes no id, or makes one twice, is refused, never the first row kept', () => {
    const d = zincFinger();
    const clash: StructureDescription = {
      ...d,
      chains: {
        label_asym_id: ['A', 'A', 'A B', 'A-1_555'],
        entity: Int32Array.of(0, 0, 0, 0),
        operator: ['1_555', '1_555', '1_555', undefined],
      },
    };
    const got = read(clash);
    expect(got.systems).toEqual([]);
    expect(got.issues.filter((i) => i.code === 'seq.io.structure.id').length).toBe(3);
    // the values of a refused chain go to no unit
    expect(got.values('b_factor')).toEqual(new Map());
  });

  it('a unit’s values are its own chain’s rows, by its key, whatever another chain says', () => {
    const d = zincFinger();
    const v = read(d).values('b_factor')!;
    expect([...v.keys()]).toEqual(['A-1_555', 'B-1_555', 'C-1_555']);
    expect(Array.from(v.get('A-1_555')!)).toEqual([Number.NaN, 20, 21, 22, 23]);
    expect(Array.from(v.get('B-1_555')!)).toEqual([30, 31, 32, 33, 34]);
    expect(read(d).values('plddt')).toBeUndefined();
    // two chains A: both refused, so B's values never land on A's unit
    const twice = { ...d, chains: { ...d.chains, label_asym_id: ['A', 'A', 'C', 'D'] } };
    expect([...read(twice).values('b_factor')!.keys()]).toEqual(['C-1_555']);
  });

  it('a position listed twice in a unit is left out (no row wins by order)', () => {
    const d = zincFinger();
    const twice: StructureDescription = {
      ...d,
      residues: { ...d.residues, position: Int32Array.of(1, 2, 3, 3, 0, 1, 2, 3, 4, 0, 0) },
      values: undefined,
    };
    const got = read(twice);
    const a = units(got.nodes.find((n) => n.key === 'complex')!.node)[0];
    expect(Array.from((a.params?.observed as { start: Int32Array }).start)).toEqual([1]);
    expect(got.issues.map((i) => i.code)).toContain('seq.io.structure.rows');
  });
});

describe('the configuration', () => {
  it('lists every limit in the Seq io table of _spec', () => {
    const spec = readFileSync(fileURLToPath(new URL('../../../_spec/README.md', import.meta.url)), 'utf8');
    const section = spec.slice(spec.indexOf('## Seq io'), spec.indexOf('## Resolve'));
    for (const key of Object.keys(DefaultStructureLimits)) expect(section).toContain(`| \`limits.${key}\``);
  });
});
