/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '@molstar/model/model/structure';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { AccessibleSurfaceArea, type ShrakeRupleyComputationProps } from '../shrake-rupley.js';

async function load1cbs() {
  const data = fs.readFileSync(
    path.resolve(import.meta.dirname, '../../../../../../../data/examples/1cbs_updated.cif'),
    'utf8',
  );
  const parsed = await CIF.parse(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

/** FNV-1a over the float32 bits: equal hashes mean bit-identical areas. */
function hash(area: ArrayLike<number>) {
  return hashFnv32a(new Uint32Array(Float32Array.from(area).buffer)).toString(16);
}

function summary(asa: AccessibleSurfaceArea) {
  let total = 0;
  for (let i = 0; i < asa.area.length; ++i) total += asa.area[i];
  return { residues: asa.area.length, total: Math.round(total), hash: hash(asa.area) };
}

describe('Shrake-Rupley accessible surface area', async () => {
  const model = await load1cbs();
  // two unit cells along a: 8 copies, units under different operators
  const cells = await StructureSymmetry.buildSymmetryRange(model, Vec3.create(0, 0, 0), Vec3.create(1, 0, 0)).run();

  const cases: [string, Structure, Partial<ShrakeRupleyComputationProps>][] = [
    ['model, defaults', model, {}],
    ['model, non-polymer occluders', model, { nonPolymer: true }],
    ['model, trace only', model, { traceOnly: true, probeSize: 4 }],
    ['model, 30 sphere points', model, { numberOfSpherePoints: 30 }],
    ['two unit cells, defaults', cells, {}],
    ['two unit cells, non-polymer occluders', cells, { nonPolymer: true }],
  ];

  // pinned from the implementation before the allocation-free computeRange (area.ts)
  const expected: Record<string, ReturnType<typeof summary>> = {
    'model, defaults': { residues: 238, total: 7782, hash: '1cffbd5d' },
    'model, non-polymer occluders': { residues: 238, total: 7561, hash: 'd1cdb060' },
    'model, trace only': { residues: 238, total: 6754, hash: '44d6afdf' },
    'model, 30 sphere points': { residues: 238, total: 7747, hash: '3eef0033' },
    'two unit cells, defaults': { residues: 1904, total: 56345, hash: 'f2c2e934' },
    'two unit cells, non-polymer occluders': { residues: 1904, total: 54474, hash: 'eb7cf8ed' },
  };

  for (const [name, structure, props] of cases) {
    it(`pins the per-residue areas: ${name}`, async () => {
      const asa = await AccessibleSurfaceArea.compute(structure, props).run();
      expect(summary(asa)).toEqual(expected[name]);
    });
  }

  it('gives the same areas when two computations interleave', async () => {
    // an update rate below zero makes every chunk notify the observer and yield (computeArea, observable.ts)
    const order: string[] = [];
    const [a, b] = await Promise.all([
      AccessibleSurfaceArea.compute(cells, {}).run(() => order.push('a'), -1),
      AccessibleSurfaceArea.compute(cells, { nonPolymer: true }).run(() => order.push('b'), -1),
    ]);
    expect(order.join('')).toMatch(/ab.*a|ba.*b/);
    expect(summary(a)).toEqual(expected['two unit cells, defaults']);
    expect(summary(b)).toEqual(expected['two unit cells, non-polymer occluders']);
  });
});

// two K+ ions (vdW 2.75 Å) 7.5 Å apart: their probe spheres overlap (2 · 1.4 + 2 · 2.75 = 8.3 Å), though they are
// farther than 2 · probe + 2 · the largest standard radius (6.54 Å) (perf-audit-2 06/B6)
const TwoIons = `data_ions
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
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_model_num
HETATM 1 K K K A 1 . 0.0 0.0 0.0 1 1
HETATM 2 K K K B 1 . 7.5 0.0 0.0 2 1
`;

describe('Shrake-Rupley lookup radius', () => {
  it('reaches the neighbours of atoms with radii above the standard ones', async () => {
    const parsed = await CIF.parse(TwoIons).run();
    if (parsed.isError) throw new Error(parsed.message);
    const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
    const structure = Structure.ofModel(trajectory.representative);
    const asa = await AccessibleSurfaceArea.compute(structure, { nonPolymer: true }).run();
    const r = 1.4 + 2.75;
    const free = 4 * Math.PI * r * r;
    expect(asa.area.length).toBe(2);
    for (const a of asa.area) {
      expect(a).toBeGreaterThan(0);
      expect(a).toBeLessThan(0.99 * free);
    }
  });
});
