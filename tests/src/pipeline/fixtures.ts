/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Large fixtures for the pipeline bench (bench.ts), made locally from 1tup.bcif (5,828 atoms, 9 units, C 1 2 1: 4
// operators per cell with the centering): its crystal operators over a block of cells, written two ways.
//
//   1tup-asm-<N>.cif    the asymmetric unit plus one assembly of N operators (pdbx_struct_oper_list): units share
//                        symmetry groups (9 groups), the capsid / large-assembly case
//   1tup-flat-<N>.bcif  the same N copies written as one model: every unit is its own group, nothing is shared (a
//                        large deposited entry, a cellPACK-style model)
//                        label_asym_id A–M repeats in every copy (non-standard; units still split per copy)
//
// usage: node scripts/pipeline-bench.mjs fixtures [--out <dir>] [--sizes 2x2x2,4x4x3,6x6x6] [--text]

import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import path from 'node:path';
import type { SymmetryOperator } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '@molstar/model/model/structure';
import { to_mmCIF } from '@molstar/model/model/structure/export/mmcif';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { ROOT, SOURCE_1TUP } from './common.js';

const args = process.argv.slice(2);
const arg = (n: string, d?: string) => {
  const i = args.indexOf(`--${n}`);
  return i >= 0 ? args[i + 1] : d;
};
const outDir = arg('out', path.join(ROOT, 'tests/build/pipeline/fixtures'))!;
const sizes = arg('sizes', '2x2x2,4x4x3,6x6x6')!.split(',');
const withText = args.includes('--text');

/** One `pdbx_struct_assembly` '1' applying every operator to every chain */
function assemblyCategories(operators: SymmetryOperator[], asymIds: string[]) {
  const lines = [
    '#',
    'loop_',
    '_pdbx_struct_assembly.id',
    '_pdbx_struct_assembly.details',
    '_pdbx_struct_assembly.method_details',
    '_pdbx_struct_assembly.oligomeric_details',
    '_pdbx_struct_assembly.oligomeric_count',
    `1 'bench assembly' ? ? ${operators.length * asymIds.length}`,
    '#',
    'loop_',
    '_pdbx_struct_assembly_gen.assembly_id',
    '_pdbx_struct_assembly_gen.oper_expression',
    '_pdbx_struct_assembly_gen.asym_id_list',
    `1 '(1-${operators.length})' ${asymIds.join(',')}`,
    '#',
    'loop_',
    '_pdbx_struct_oper_list.id',
    '_pdbx_struct_oper_list.type',
    '_pdbx_struct_oper_list.name',
    '_pdbx_struct_oper_list.matrix[1][1]',
    '_pdbx_struct_oper_list.matrix[1][2]',
    '_pdbx_struct_oper_list.matrix[1][3]',
    '_pdbx_struct_oper_list.vector[1]',
    '_pdbx_struct_oper_list.matrix[2][1]',
    '_pdbx_struct_oper_list.matrix[2][2]',
    '_pdbx_struct_oper_list.matrix[2][3]',
    '_pdbx_struct_oper_list.vector[2]',
    '_pdbx_struct_oper_list.matrix[3][1]',
    '_pdbx_struct_oper_list.matrix[3][2]',
    '_pdbx_struct_oper_list.matrix[3][3]',
    '_pdbx_struct_oper_list.vector[3]',
  ];
  const f = (x: number) => +x.toFixed(6);
  operators.forEach((op, i) => {
    const m = op.matrix; // column-major
    const row = (r: number) => `${f(m[r])} ${f(m[r + 4])} ${f(m[r + 8])} ${f(m[r + 12])}`;
    lines.push(`${i + 1} 'crystal symmetry operation' ${op.name} ${row(0)} ${row(1)} ${row(2)}`);
  });
  lines.push('#', '');
  return lines.join('\n');
}

const parsed = await CIF.parseBinary(new Uint8Array(readFileSync(path.join(ROOT, SOURCE_1TUP)))).run();
if (parsed.isError) throw new Error(parsed.message);
const model = (await trajectoryFromMmCIF(parsed.result.blocks[0]).run()).representative;
const asym = Structure.ofModel(model);
const asymIds = model.atomicHierarchy.chains.label_asym_id.toArray() as string[];
mkdirSync(outDir, { recursive: true });

for (const size of sizes) {
  const [i, j, k] = size.split('x').map((x) => +x);
  const sym = await StructureSymmetry.buildSymmetryRange(
    asym,
    Vec3.create(0, 0, 0),
    Vec3.create(i - 1, j - 1, k - 1),
  ).run();
  const operators = new Map<string, SymmetryOperator>();
  for (const u of sym.units) operators.set(u.conformation.operator.name, u.conformation.operator);
  const n = operators.size;

  const asymCif = to_mmCIF('1TUP_ASM', asym, false, {
    skipCategoryNames: new Set(['pdbx_struct_assembly', 'pdbx_struct_assembly_gen', 'pdbx_struct_oper_list']),
  }) as string;
  writeFileSync(path.join(outDir, `1tup-asm-${n}.cif`), asymCif + assemblyCategories([...operators.values()], asymIds));

  writeFileSync(path.join(outDir, `1tup-flat-${n}.bcif`), to_mmCIF('1TUP_FLAT', sym, true) as Uint8Array);
  if (withText) writeFileSync(path.join(outDir, `1tup-flat-${n}.cif`), to_mmCIF('1TUP_FLAT', sym, false) as string);
  console.log(
    `${size}: ${n} operators, ${sym.elementCount} atoms, ${sym.units.length} units -> 1tup-asm-${n}.cif, 1tup-flat-${n}.bcif`,
  );
}
