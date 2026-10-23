/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { type ElementIndex, Structure, StructureElement } from '@molstar/model/model/structure';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { SortedArray } from '@molstar/core/data/int';
import { AccessibleSurfaceArea } from '../accessible-surface-area/shrake-rupley.js';
import { AccessibleSurfaceAreaProvider } from '../accessible-surface-area.js';

const Residue = [
  ['N', 0.0, 1.4, 0.0],
  ['CA', 0.0, 0.0, 0.0],
  ['C', 1.5, -0.5, 0.0],
  ['O', 1.8, -1.7, 0.0],
  ['CB', -0.8, -0.6, 1.2],
] as const;

/** A 3-residue peptide; assembly 1 is 3 copies: as is, 4.5 Å away in z (in contact) and 100 Å away. */
function mmcif() {
  const atoms: string[] = [];
  let id = 1;
  for (let r = 0; r < 3; r++) {
    for (const [name, x, y, z] of Residue) {
      atoms.push(
        `ATOM ${id++} ${name[0]} ${name} ALA A 1 ${r + 1} ${x + 3.8 * r} ${y} ${z} 1 10 ${r + 1} ALA A ${name} 1`,
      );
    }
  }
  const opers = [0, 4.5, 100]
    .map((dz, i) => `${i + 1} 'point symmetry operation' 1 0 0 0 0 1 0 0 0 0 1 ${dz}`)
    .join('\n');
  return `data_T
_pdbx_struct_assembly.id 1
_pdbx_struct_assembly.details author_defined_assembly
_pdbx_struct_assembly_gen.assembly_id 1
_pdbx_struct_assembly_gen.oper_expression (1,2,3)
_pdbx_struct_assembly_gen.asym_id_list A
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
${opers}
loop_
${'group_PDB id type_symbol label_atom_id label_comp_id label_asym_id label_entity_id label_seq_id Cartn_x Cartn_y Cartn_z occupancy B_iso_or_equiv auth_seq_id auth_comp_id auth_asym_id auth_atom_id pdbx_PDB_model_num'
  .split(' ')
  .map((c) => `_atom_site.${c}`)
  .join('\n')}
${atoms.join('\n')}
`;
}

async function assembly() {
  const parsed = await CIF.parseText(mmcif()).run();
  if (parsed.isError) throw new Error(parsed.message);
  const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
  return await StructureSymmetry.buildAssembly(Structure.ofModel(model), '1').run();
}

/** The value of each element of each unit, in the order the computation stores them (by serial index of `root`). */
function expected(root: Structure, asa: AccessibleSurfaceArea) {
  const values = new Map<number, number[]>();
  let m = 0;
  for (const unit of root.units) {
    const xs: number[] = [];
    for (let j = 0; j < unit.elements.length; j++) {
      const rSI = asa.serialResidueIndex[m + j];
      xs.push(rSI === -1 ? -1 : asa.area[rSI]);
    }
    values.set(unit.id, xs);
    m += unit.elements.length;
  }
  return values;
}

function values(structure: Structure, asa: AccessibleSurfaceArea) {
  const l = StructureElement.Location.create(structure);
  const ret = new Map<number, number[]>();
  for (const unit of structure.units) {
    l.unit = unit;
    const xs: number[] = [];
    for (let j = 0; j < unit.elements.length; j++) {
      l.element = unit.elements[j];
      xs.push(AccessibleSurfaceArea.getValue(l, asa));
    }
    ret.set(unit.id, xs);
  }
  return ret;
}

describe('AccessibleSurfaceArea.getValue', () => {
  it('reads each copy of an assembly from its own unit', async () => {
    const structure = await assembly();
    expect(structure.units.length).toBe(3);
    const asa = await AccessibleSurfaceArea.compute(structure).run();

    const exp = expected(structure, asa);
    expect(values(structure, asa)).toEqual(exp);

    // the copies in contact bury some of each other, the far one does not
    const [a, b, c] = structure.units.map((u) => exp.get(u.id)!);
    expect(a).not.toEqual(c);
    expect(b).not.toEqual(c);
    for (const xs of [a, b, c]) for (const x of xs) expect(x).toBeGreaterThan(0);
  });

  it('reads a substructure from its root', async () => {
    const root = await assembly();
    const asa = await AccessibleSurfaceArea.compute(root).run();
    const exp = expected(root, asa);

    // the second residue of the second and third copies
    const units = root.units.slice(1).map((u) => u.getChild(SortedArray.ofRange(5 as ElementIndex, 9 as ElementIndex)));
    const sub = Structure.create(units, { parent: root });
    expect(sub.root).toBe(root);

    const got = values(sub, asa);
    for (const unit of sub.units) {
      expect(got.get(unit.id)).toEqual(exp.get(unit.id)!.slice(5, 10));
    }
  });
});

describe('AccessibleSurfaceArea on coarse units', () => {
  it('gives coarse elements no value instead of throwing, and does not apply to coarse-only structures', async () => {
    const examples = path.resolve(import.meta.dirname, '../../../../../../data/examples');
    const parsed = await CIF.parseText(fs.readFileSync(path.join(examples, 'ihm-coarse-8zzh.cif'), 'utf8')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
    const structure = Structure.ofModel(model);
    expect(structure.units.length).toBeGreaterThan(0);

    const asa = await AccessibleSurfaceArea.compute(structure).run();
    expect(asa.area.length).toBe(0);
    const l = StructureElement.Location.create(structure, structure.units[0], structure.units[0].elements[0]);
    expect(AccessibleSurfaceArea.getFlag(l, asa)).toBe(AccessibleSurfaceArea.Flags.NA);

    expect(AccessibleSurfaceAreaProvider.isApplicable(structure)).toBe(false);
    expect(AccessibleSurfaceAreaProvider.isApplicable(await assembly())).toBe(true);
  });
});
