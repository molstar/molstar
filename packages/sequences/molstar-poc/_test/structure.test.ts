/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `molstar:structure` over the committed fixtures (seqspec:libraries#the-mol-adapter): units and numberings of each,
// SIFTS, the options, and a structure's units read through the load path (services, formats, a by_numbering mapping
// onto P04637).

import { compile } from '@molstar/seq-core/format';
import { KeyKind, type NumberingInput, type SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { complexInputs, readStructure, type StructureOptions, structureOf } from '../src/index.js';
import { bytesOf, loaders, P53, p53Doc, registries } from './fixtures.js';

async function inputs(file: string, options: Partial<StructureOptions> = {}) {
  const src = await readStructure(bytesOf(file), 'bcif');
  const structure = await structureOf(src, 0, options.assembly ?? null);
  const issues: string[] = [];
  const r = complexInputs([{ structure, model: null }], src.frame, options, (code) => issues.push(code));
  return { ...r, issues };
}
const labelsOf = (n: NumberingInput | undefined) =>
  n?.method === 'labels' ? n.labels : n?.method === 'names' ? n.names : null;

describe('molstar:structure', () => {
  it('1TUP: DNA, p53 and zinc as molecules; chains as units with observed residues, auth labels and mvs fields', async () => {
    const r = await inputs('1tup.bcif');
    expect(r.issues).toEqual([]);
    expect(
      r.molecules.map((m) => [m.key, m.kind ?? 'polymer', m.alphabet ?? null, m.residues?.length ?? null]),
    ).toEqual([
      ['entity:1', 'polymer', 'dna', 21],
      ['entity:2', 'polymer', 'dna', 21],
      ['entity:3', 'polymer', 'protein', 219],
      ['entity:4', 'ligand', null, null],
    ]);
    const p53 = r.molecules[2];
    // The construct is P04637 94–312, letter for letter.
    expect(p53.residues).toBe(P53.sequence.slice(93, 312));
    expect(p53.numberings).toEqual([
      { name: 'label', method: 'offset', first: 1 },
      { name: 'uniprot:P04637', method: 'segments', segments: [[0, 94, 219]] },
    ]);
    expect(p53.identifiers).toEqual([{ namespace: 'uniprot', value: 'P04637' }]);
    expect(r.units.map((u) => [u.key, u.molecule])).toEqual([
      ['A', 'entity:1'],
      ['B', 'entity:2'],
      ['C', 'entity:3'],
      ['D', 'entity:3'],
      ['E', 'entity:3'],
      ['F', 'entity:4'],
      ['G', 'entity:4'],
      ['H', 'entity:4'],
    ]);
    const [c, d] = [r.units[2], r.units[3]];
    expect([c.observed, d.observed]).toEqual([[[0, 196]], [[2, 196]]]);
    expect(c.custom).toEqual({ mvs: { label_asym_id: 'C', auth_asym_id: 'A', instance_id: '1_555' } });
    expect(c.origin).toBe('experimental');
    // auth = UniProt here; unobserved residues (196–218) keep pdbx_poly_seq_scheme's numbers.
    const auth = labelsOf(c.numberings?.[0]);
    expect([auth?.[0], auth?.[81], auth?.[218], c.numberings?.[0].isDefault]).toEqual([94, 175, 312, true]);
    expect(r.units[5].numberings).toEqual([{ name: 'atom', method: 'names', names: ['ZN'], isDefault: true }]);
  });

  it('2OCJ and 2XWR: SIFTS segments onto P04637; the AlphaFold model is predicted, P04637 whole', async () => {
    const ocj = await inputs('2ocj.bcif');
    expect(ocj.molecules[0].numberings?.[1]).toEqual({
      name: 'uniprot:P04637',
      method: 'segments',
      segments: [[0, 94, 219]],
    });
    expect(ocj.units.filter((u) => u.molecule === 'entity:1').map((u) => u.key)).toEqual(['A', 'B', 'C', 'D']);
    const xwr = await inputs('2xwr.bcif');
    expect(xwr.molecules[0].numberings?.[1]).toEqual({
      name: 'uniprot:P04637',
      method: 'segments',
      segments: [[0, 89, 205]],
    });
    expect(xwr.molecules[0].residues).toBe(P53.sequence.slice(88, 293));
    const af = await inputs('AF-P04637-F1-model_v6.bcif');
    expect(af.molecules.map((m) => [m.key, m.residues])).toEqual([['entity:1', P53.sequence]]);
    expect(af.units.map((u) => [u.key, u.origin, u.observed ?? null])).toEqual([['A', 'predicted', null]]);
  });

  it("options: SIFTS isoforms or none, chains, include, an assembly's copies, unit-name collisions", async () => {
    const all = await inputs('1tup.bcif', { sifts: 'all' });
    expect(all.molecules[2].numberings?.map((n) => n.name)).toEqual([
      'label',
      ...['P04637', ...[2, 3, 4, 5, 6, 7, 8, 9].map((k) => `P04637-${k}`)].map((a) => `uniprot:${a}`),
    ]);
    expect(all.molecules[2].numberings?.find((n) => n.name === 'uniprot:P04637-7')).toMatchObject({
      segments: [[39, 1, 180]],
    });
    expect((await inputs('1tup.bcif', { sifts: 'none' })).molecules[2].numberings?.map((n) => n.name)).toEqual([
      'label',
    ]);
    const some = await inputs('1tup.bcif', { chains: ['C', 'F'], include: ['polymer'] });
    expect([some.units.map((u) => u.key), some.molecules.map((m) => m.key)]).toEqual([['C'], ['entity:3']]);
    // Auth ids collide where an ion shares its chain's; a template tells them apart.
    expect((await inputs('1tup.bcif', { unit_name: 'auth_asym_id' })).issues).toEqual(['molstar.unit_name.collision']);
    const named = await inputs('1tup.bcif', { unit_name: '{auth_asym_id}-{entity_id}' });
    expect(named.units.map((u) => u.key)).toEqual(['E-1', 'F-2', 'A-3', 'B-3', 'C-3', 'A-4', 'B-4', 'C-4']);
    // 1C26's tetramer: chain A four times, each copy named by its operator off the identity.
    const tetramer = await inputs('1c26.bcif', { assembly: '1' });
    expect(tetramer.units.length).toBe(4);
    expect(new Set(tetramer.units.map((u) => u.key)).size).toBe(4);
    expect(tetramer.units.every((u) => u.key === 'A' || u.key.startsWith('A@'))).toBe(true);
    expect(tetramer.molecules.length).toBe(1);
  });
});

const unitId = (set: SeqSet, ref: string) => {
  const t = { kind: 0, id: 0, row: -1 };
  return set.byName(ref, t) && t.kind === KeyKind.unit ? t.id : -1;
};

describe('through the load path', () => {
  it("pdbe and uniprot sources, the complex's units, and a by_numbering mapping from P04637 onto 1TUP's p53", async () => {
    const { set, issues } = await compile(p53Doc([{ ref: '1tup', pdb: '1TUP', entity: '3' }]), registries(), loaders());
    expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
    expect(unitId(set, '1tup/C')).toBeGreaterThanOrEqual(0);
    const t = { kind: 0, id: 0, row: -1 };
    set.byName('1tup-map', t);
    const m = set.mappings[t.id];
    // R175 (index 174) is the construct's 82nd residue (index 81); outside 94–312 nothing.
    expect([m?.forward[174], m?.forward[93], m?.forward[92], m?.forward[312]]).toEqual([81, 0, -1, -1]);
  });
});
