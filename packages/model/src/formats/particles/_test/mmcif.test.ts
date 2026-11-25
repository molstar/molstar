/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { CIF, CifBlock } from '@molstar/io/reader/cif';
import { toDatabase } from '@molstar/io/reader/cif/schema';
import { mmCIF_Schema } from '@molstar/io/reader/cif/schema/mmcif';
import { getMatrices, expandOperators, parseOperatorList } from '../../structure/property/assembly.js';
import {
  composeOperatorCombination,
  expandOperatorCombinations,
  formatOperatorCombination,
  getParticleOperator,
  getParticleOperatorId,
  getParticleOperatorIndex,
  parseOperatorCombinations,
  readParticleOperators,
} from '../operators.js';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { type ParticleList, ParticleTarget, getParticleTransforms } from '@molstar/model/model/particles/particle-list';
import { createParticleListFromMmcifAssembly } from '../mmcif.js';
import { RuntimeContext, Task } from '@molstar/core/task';
import { ElementSymbol } from '@molstar/model/model/structure/model/types';
import { VdwRadius } from '@molstar/model/model/structure/model/properties/atomic/measures';

const rows = [
  // id, matrix (row-major 3x3), vector
  ['1', [1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0]],
  ['2', [0, -1, 0, 1, 0, 0, 0, 0, 1], [10, 20, 30]],
  ['P', [0.5, 0.5, -0.7071, -0.5, 0.8536, 0.1464, 0.7071, 0.1464, 0.6929], [-1.5, 2.25, 3]],
] as const;

function makeCif(naming: 'brackets' | 'brackets0' | 'underscore', ids?: readonly string[]) {
  const names: string[] = ['_pdbx_struct_oper_list.id'];
  const off = naming === 'brackets0' ? 0 : 1;
  for (let i = 0; i < 3; ++i) {
    for (let j = 0; j < 3; ++j) {
      names.push(
        naming === 'underscore'
          ? `_pdbx_struct_oper_list.matrix_${i + 1}${j + 1}`
          : `_pdbx_struct_oper_list.matrix[${i + off}][${j + off}]`,
      );
    }
  }
  for (let i = 0; i < 3; ++i) {
    names.push(
      naming === 'underscore' ? `_pdbx_struct_oper_list.vector_${i + 1}` : `_pdbx_struct_oper_list.vector[${i + off}]`,
    );
  }

  const lines = ['data_test', 'loop_', ...names];
  for (let i = 0; i < rows.length; ++i) {
    const [id, matrix, vector] = rows[i];
    lines.push([ids ? ids[i] : id, ...matrix, ...vector].join(' '));
  }
  return lines.join('\n') + '\n';
}

async function parseBlock(naming: 'brackets' | 'brackets0' | 'underscore', ids?: readonly string[]): Promise<CifBlock> {
  const parsed = await CIF.parseText(makeCif(naming, ids)).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result.blocks[0];
}

const cellpackFiberCif = `data_cellpack_test
loop_
_entity.id
_entity.pdbx_description
1 root.test.fibers.F
2 root.test.proteins.P
3 root.test.fiber.N

loop_
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
F 1 0 0 0
P 2 0 0 0
N 3 0 0 0

loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(3,1,2)' F
1 '(4-5)' F
1 '(6)' F
1 '(1-2)' P
1 '(1-2)' N

loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.vector[3]
1 1 0 0 0 1 0 0 0 1 10 0 0
2 1 0 0 0 1 0 0 0 1 20 0 0
3 1 0 0 0 1 0 0 0 1 30 0 0
4 1 0 0 0 1 0 0 0 1 40 0 0
5 1 0 0 0 1 0 0 0 1 50 0 0
6 1 0 0 0 1 0 0 0 1 60 0 0
`;

async function parseCellpackFiberCif() {
  const parsed = await CIF.parseText(cellpackFiberCif).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result;
}

describe('mmcif particle operators', () => {
  for (const naming of ['brackets', 'brackets0', 'underscore'] as const) {
    it(`matches getMatrices (${naming})`, async () => {
      const block = await parseBlock(naming);
      const db = toDatabase(mmCIF_Schema, block);
      const expected = getMatrices(db.pdbx_struct_oper_list);

      const ops = readParticleOperators(block);
      expect(ops.count).toBe(rows.length);
      expect(rows.map((_, i) => getParticleOperatorId(ops, i))).toEqual(rows.map((r) => r[0]));

      const m = Mat4();
      for (let i = 0; i < ops.count; ++i) {
        getParticleOperator(ops, i, m);
        const ref = expected.get(getParticleOperatorId(ops, i))!;
        for (let k = 0; k < 16; ++k) {
          expect(m[k]).toBeCloseTo(ref[k], 10);
        }
      }
    });
  }

  it('skips the id string table for serial 1-based ids', async () => {
    const serial = await parseBlock('brackets', ['1', '2', '3']);
    const ops = readParticleOperators(serial);
    expect(ops.ids).toBeUndefined();
    expect(ops.index).toBeUndefined();
    expect(getParticleOperatorId(ops, 2)).toBe('3');
    expect(expandOperatorCombinations(parseOperatorList('(3,1)'), ops).indices).toEqual(new Int32Array([2, 0]));
    expect(() => expandOperatorCombinations(parseOperatorList('(4)'), ops)).toThrow(/Operator '4' not found/);

    const named = await parseBlock('brackets');
    expect(readParticleOperators(named).ids).toEqual(['1', '2', 'P']);
  });

  it('stores integer ids with gaps without a string table', async () => {
    const block = await parseBlock('brackets', ['1', '3', '10']);
    const ops = readParticleOperators(block);
    expect(ops.ids).toEqual(new Int32Array([1, 3, 10]));
    expect(getParticleOperatorId(ops, 2)).toBe('10');
    expect(getParticleOperatorIndex(ops, '3')).toBe(1);
    expect(getParticleOperatorIndex(ops, '2')).toBe(-1);
    expect(getParticleOperatorIndex(ops, '03')).toBe(-1);
    expect(expandOperatorCombinations(parseOperatorList('(10,1)'), ops).indices).toEqual(new Int32Array([2, 0]));
  });

  it('expands combinations in the same order as expandOperators', async () => {
    const block = await parseBlock('brackets');
    const ops = readParticleOperators(block);

    for (const expression of ['(1,2,P)', '(1,2)(2,P)', '(P)']) {
      const list = parseOperatorList(expression);
      const expected = expandOperators(list);
      const combinations = expandOperatorCombinations(list, ops);

      expect(combinations.count).toBe(expected.length);
      expect(combinations.stride).toBe(list.length);

      for (let i = 0; i < combinations.count; ++i) {
        expect(formatOperatorCombination(ops.ids, combinations, i)).toBe(expected[i].join('×'));
      }
    }
  });

  it('parses expressions to the combinations of the operator lists', async () => {
    const serial = readParticleOperators(await parseBlock('brackets', ['1', '2', '3']));
    const gaps = readParticleOperators(await parseBlock('brackets', ['1', '3', '10']));
    const named = readParticleOperators(await parseBlock('brackets'));
    for (const [ops, expressions] of [
      [serial, ['(1-3)', '(1,2,3)(2-3)', '(3)', '1-2', '(3-1)']],
      [gaps, ['(1,3,10)', '(10)(1,3)', '(3-3)']],
      [named, ['(1,2,P)', '(1-2)(P)', '(P)']],
    ] as const) {
      for (const expression of expressions) {
        expect(parseOperatorCombinations(expression, ops)).toEqual(
          expandOperatorCombinations(parseOperatorList(expression), ops),
        );
      }
    }
    expect(() => parseOperatorCombinations('(1-4)', serial)).toThrow(/Operator '4' not found/);
    expect(() => parseOperatorCombinations('(2-3)', gaps)).toThrow(/Operator '2' not found/);
    expect(() => parseOperatorCombinations('(X)', named)).toThrow(/Operator 'X' not found/);
  });

  it('composes combinations like sequential Mat4.mul', async () => {
    const block = await parseBlock('brackets');
    const db = toDatabase(mmCIF_Schema, block);
    const matrices = getMatrices(db.pdbx_struct_oper_list);
    const ops = readParticleOperators(block);

    const list = parseOperatorList('(1,2,P)(2,P)');
    const expected = expandOperators(list);
    const combinations = expandOperatorCombinations(list, ops);

    const actual = Mat4();
    const ref = Mat4();
    for (let i = 0; i < combinations.count; ++i) {
      composeOperatorCombination(ops, combinations, i, actual);

      Mat4.setIdentity(ref);
      for (const id of expected[i]) Mat4.mul(ref, ref, matrices.get(id)!);

      for (let k = 0; k < 16; ++k) {
        expect(actual[k]).toBeCloseTo(ref[k], 10);
      }
    }
  });

  it('throws for unknown operator ids', async () => {
    const block = await parseBlock('brackets');
    const ops = readParticleOperators(block);
    expect(() => expandOperatorCombinations(parseOperatorList('(9)'), ops)).toThrow(/Operator '9' not found/);
  });

  it('returns an empty table when the category is missing', async () => {
    const parsed = await CIF.parseText('data_test\n_entry.id test\n').run();
    if (parsed.isError) throw new Error(parsed.message);
    const ops = readParticleOperators(parsed.result.blocks[0]);
    expect(ops.count).toBe(0);
    expect(ops.columns.length).toBe(0);
  });
});

describe('CellPack mmCIF fibers', () => {
  it('uses separate assembly-gen rows as ordered fibers', async () => {
    const cif = await parseCellpackFiberCif();
    const list = await createParticleListFromMmcifAssembly(RuntimeContext.Synchronous, cif, {
      assemblyId: '1',
      variant: 'cellpack',
    });

    expect(list.count).toBe(10);
    expect(list.fibers).toBeDefined();
    expect(list.fibers!.count).toBe(2);
    expect(Array.from(list.fibers!.offsets)).toEqual([0, 3, 5]);
    expect(Array.from(list.fibers!.indices)).toEqual([0, 1, 2, 3, 4]);
    expect([list.coordinates[0], list.coordinates[3], list.coordinates[6]]).toEqual([30, 10, 20]);
  });

  it('respects asym filtering and omits singleton fibers', async () => {
    const cif = await parseCellpackFiberCif();
    const list = await createParticleListFromMmcifAssembly(RuntimeContext.Synchronous, cif, {
      assemblyId: '1',
      asymIds: ['F'],
      variant: 'cellpack',
    });

    expect(list.count).toBe(6);
    expect(list.fibers!.count).toBe(2);
    expect(Array.from(list.fibers!.offsets)).toEqual([0, 3, 5]);
  });

  it('does not infer fibers for ordinary entities or the standard variant', async () => {
    const cif = await parseCellpackFiberCif();
    const ordinary = await createParticleListFromMmcifAssembly(RuntimeContext.Synchronous, cif, {
      assemblyId: '1',
      asymIds: ['P', 'N'],
      variant: 'cellpack',
    });
    const standard = await createParticleListFromMmcifAssembly(RuntimeContext.Synchronous, cif, {
      assemblyId: '1',
      asymIds: ['F'],
      variant: 'standard',
    });

    expect(ordinary.fibers).toBeUndefined();
    expect(standard.fibers).toBeUndefined();
  });

  it('reports progress while calculating bounds and building particles', async () => {
    const cif = await parseCellpackFiberCif();
    const updates: string[] = [];
    const task = Task.create('Create particle list', (ctx) =>
      createParticleListFromMmcifAssembly(ctx, cif, { assemblyId: '1', variant: 'cellpack' }),
    );
    const list = await task.run((progress) => updates.push(progress.root.progress.message), 0);

    expect(list.count).toBe(10);
    expect(updates).toContain('Reading particle targets');
    expect(updates).toContain('Building particles');
  });
});

// one chain whose centroid (11/3, 0, 0) is not its box's center (5, 0, 0), under a rotation and a translation
const placementCif = `data_placement_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA ALA A 1 1 0 0 0 1
ATOM 2 C CA ALA A 1 2 1 0 0 1
ATOM 3 C CA ALA A 1 3 10 0 0 1

loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(1-2)' A

loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.vector[3]
1 1 0 0 0 1 0 0 0 1 0 0 0
2 0 -1 0 1 0 0 0 0 1 10 20 30
`;

// PetWorld: a model a target; model 1's chains out of a model's order (A of entity 1 interleaved, entity 2 between)
const petworldCif = (asym: 'label_asym_id' | 'auth_asym_id') => `data_petworld_chains
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.${asym}
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA ALA A 1 1 0 0 0 1
ATOM 2 C CA ALA B 2 1 1 0 0 1
ATOM 3 C CA ALA A 1 2 2 0 0 1
ATOM 4 C CA ALA C 1 1 3 0 0 1
ATOM 5 C CA ALA A 1 1 5 0 0 2
ATOM 6 C CA ALA B 1 1 6 0 0 2

loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
_pdbx_struct_assembly_gen.PDB_model_num
1 '(1)' A 1
1 '(1)' A 2

loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.vector[3]
1 1 0 0 0 1 0 0 0 1 0 0 0
`;

async function parsePetworld(text: string) {
  const parsed = await CIF.parseText(text).run();
  if (parsed.isError) throw new Error(parsed.message);
  return createParticleListFromMmcifAssembly(RuntimeContext.Synchronous, parsed.result, { assemblyId: '1' });
}

function sitesOf(list: ParticleList, id: number) {
  const target = list.targetMapping!.get(id)!;
  if (target.kind !== 'sites') throw new Error('expected a sites target');
  return target.sites;
}

describe('mmCIF particle targets', () => {
  it('place the atoms where the operators put them', async () => {
    const parsed = await CIF.parseText(placementCif).run();
    if (parsed.isError) throw new Error(parsed.message);
    const list = await createParticleListFromMmcifAssembly(RuntimeContext.Synchronous, parsed.result, {
      assemblyId: '1',
      variant: 'standard',
    });
    const target = list.targetMapping!.get(0)!;
    const center = ParticleTarget.getCenter(target)!;
    expect(center).toBe(target.kind === 'sites' ? target.center : undefined);
    for (let k = 0; k < 3; ++k) expect(center[k]).toBeCloseTo([11 / 3, 0, 0][k], 5);

    const ops = readParticleOperators(parsed.result.blocks[0]);
    const transforms = getParticleTransforms(list);
    const op = Mat4();
    const m = Mat4();
    const placed = Vec3();
    const expected = Vec3();
    for (let p = 0; p < list.count; ++p) {
      getParticleOperator(ops, p, op);
      Mat4.fromArray(m, transforms, p * 16);
      for (const x of [0, 1, 10]) {
        const atom = Vec3.create(x, 0, 0);
        Vec3.transformMat4(placed, Vec3.sub(Vec3(), atom, center), m);
        Vec3.transformMat4(expected, atom, op);
        for (let k = 0; k < 3; ++k) expect(placed[k]).toBeCloseTo(expected[k], 4);
      }
    }
  });

  it('are the sites of the chains: positions, van der Waals radii, centered by their centroid', async () => {
    const parsed = await CIF.parseText(placementCif).run();
    if (parsed.isError) throw new Error(parsed.message);
    const list = await createParticleListFromMmcifAssembly(RuntimeContext.Synchronous, parsed.result, {
      assemblyId: '1',
      variant: 'standard',
    });
    const target = list.targetMapping!.get(0)!;
    if (target.kind !== 'sites') throw new Error('expected a sites target');
    const { sites } = target;
    expect(sites.count).toBe(3);
    expect(Array.from(sites.positions)).toEqual([0, 0, 0, 1, 0, 0, 10, 0, 0]);
    // carbon's, as a model's physical size (`VdwRadiusOfNumber` of its atomic number)
    expect(Array.from(sites.radii)).toEqual(new Array(3).fill(Math.fround(VdwRadius(ElementSymbol('C')))));
    // without the given center: the sites' centroid
    const center = ParticleTarget.getCenter({ kind: 'sites', sites })!;
    for (let k = 0; k < 3; ++k) expect(center[k]).toBeCloseTo([11 / 3, 0, 0][k], 5);
  });

  it("are a chain's sites as one run, labelled by its label_asym_id", async () => {
    const parsed = await CIF.parseText(placementCif).run();
    if (parsed.isError) throw new Error(parsed.message);
    const list = await createParticleListFromMmcifAssembly(RuntimeContext.Synchronous, parsed.result, {
      assemblyId: '1',
      variant: 'standard',
    });
    const target = list.targetMapping!.get(0)!;
    if (target.kind !== 'sites') throw new Error('expected a sites target');
    expect(Array.from(target.sites.chains!.offsets)).toEqual([0, 3]);
    expect(target.sites.chains!.labels).toEqual(['A']);
  });

  it("are a model's sites by chain in a model's order: reordered only where out of it", async () => {
    const list = await parsePetworld(petworldCif('label_asym_id'));
    // model 1: A (entity 1), B (entity 2), A, C (entity 1): by entity, then chain, each as it first appears
    const first = sitesOf(list, 0);
    expect(Array.from(first.positions.filter((_, i) => i % 3 === 0))).toEqual([0, 2, 3, 1]);
    expect(Array.from(first.chains!.offsets)).toEqual([0, 2, 3, 4]);
    expect(first.chains!.labels).toEqual(['A', 'C', 'B']);
    // model 2: in order already, its rows' order kept
    const second = sitesOf(list, 1);
    expect(Array.from(second.positions.filter((_, i) => i % 3 === 0))).toEqual([5, 6]);
    expect(Array.from(second.chains!.offsets)).toEqual([0, 1, 2]);
    expect(second.chains!.labels).toEqual(['A', 'B']);
  });

  it("are a model's chains by auth_asym_id without label_asym_id", async () => {
    const list = await parsePetworld(petworldCif('auth_asym_id'));
    expect(sitesOf(list, 0).chains!.labels).toEqual(['A', 'C', 'B']);
  });

  it('are none without resolving targets', async () => {
    const parsed = await CIF.parseText(placementCif).run();
    if (parsed.isError) throw new Error(parsed.message);
    const list = await createParticleListFromMmcifAssembly(RuntimeContext.Synchronous, parsed.result, {
      assemblyId: '1',
      variant: 'standard',
      resolveTargets: false,
    });
    expect(list.targetMapping).toBeUndefined();
  });
});
