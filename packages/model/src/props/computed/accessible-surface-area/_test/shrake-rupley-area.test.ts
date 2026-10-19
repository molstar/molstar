/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { SymmetryOperator } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { AssetManager } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, type Unit } from '@molstar/model/model/structure';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { AccessibleSurfaceAreaProvider } from '../../accessible-surface-area.js';
import {
  AccessibleSurfaceArea,
  ShrakeRupleyComputationParams,
  type ShrakeRupleyComputationProps,
} from '../shrake-rupley.js';

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
    ['two unit cells, every unit', cells, { symmetry: 'all' }],
    ['two unit cells, non-polymer occluders', cells, { nonPolymer: true }],
  ];

  // pinned from the implementation before the allocation-free computeRange (area.ts)
  const expected: Record<string, ReturnType<typeof summary>> = {
    'model, defaults': { residues: 238, total: 7782, hash: '1cffbd5d' },
    'model, non-polymer occluders': { residues: 238, total: 7561, hash: 'd1cdb060' },
    'model, trace only': { residues: 238, total: 6754, hash: '44d6afdf' },
    'model, 30 sphere points': { residues: 238, total: 7747, hash: '3eef0033' },
    'two unit cells, defaults': { residues: 1904, total: 56345, hash: 'f2c2e934' },
    'two unit cells, every unit': { residues: 1904, total: 56345, hash: 'f2c2e934' },
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

/**
 * Copies of `structure` rotated by `2π k / n` for each k of `steps`, about an axis parallel to z `offset` Å from its
 * centre along x. All n steps make C_n, a closed point group: every copy has the same surroundings, rotated.
 */
function cyclicAssembly(structure: Structure, n: number, offset: number, steps: number[]) {
  const axis = Vec3.add(Vec3(), structure.boundary.sphere.center, Vec3.create(offset, 0, 0));
  const builder = Structure.Builder({ label: `C${n}` });
  for (const k of steps) {
    const m = Mat4.fromTranslation(Mat4(), axis);
    Mat4.mul(m, m, Mat4.fromRotation(Mat4(), (2 * Math.PI * k) / n, Vec3.unitZ));
    Mat4.mul(m, m, Mat4.fromTranslation(Mat4(), Vec3.negate(Vec3(), axis)));
    const op = SymmetryOperator.create(`C${n}-${k}`, m);
    for (const unit of structure.units) builder.addWithOperator(unit, op);
  }
  return builder.getStructure();
}

/** The area of each element's residue, per unit id (-1 for elements without one), by serial index of `structure` */
function areasByUnit(structure: Structure, asa: AccessibleSurfaceArea) {
  const ret = new Map<number, number[]>();
  let m = 0;
  for (const unit of structure.units) {
    const xs: number[] = [];
    for (let j = 0; j < unit.elements.length; ++j) {
      const rSI = asa.serialResidueIndex[m + j];
      xs.push(rSI === -1 ? -1 : asa.area[rSI]);
    }
    ret.set(unit.id, xs);
    m += unit.elements.length;
  }
  return ret;
}

/** The copy (of `copies`, `unitsPerCopy` units each) whose first unit has the same areas under 'per-group' as 'all' */
function representativeCopy(
  structure: Structure,
  all: Map<number, number[]>,
  perGroup: Map<number, number[]>,
  unitsPerCopy: number,
  copies: number,
) {
  const exact: number[] = [];
  for (let k = 0; k < copies; ++k) {
    const id = structure.units[k * unitsPerCopy].id;
    if (maxDifference(all.get(id)!, perGroup.get(id)!) === 0) exact.push(k);
  }
  expect(exact.length).toBe(1);
  return exact[0];
}

function maxDifference(a: number[], b: number[]) {
  let max = 0;
  for (let i = 0; i < a.length; ++i) max = Math.max(max, Math.abs(a[i] - b[i]));
  return max;
}

describe('Shrake-Rupley accessible surface area, once per symmetry group', async () => {
  const model = await load1cbs();
  const cells = await StructureSymmetry.buildSymmetryRange(model, Vec3.create(0, 0, 0), Vec3.create(1, 0, 0)).run();

  it("defaults to 'all'", () => {
    expect(PD.getDefaultValues(ShrakeRupleyComputationParams).symmetry).toBe('all');
  });

  it("is the same as 'all' on a structure without copies", async () => {
    const all = await AccessibleSurfaceArea.compute(model, { symmetry: 'all' }).run();
    const perGroup = await AccessibleSurfaceArea.compute(model, { symmetry: 'per-group' }).run();
    expect(summary(perGroup)).toEqual(summary(all));
  });

  it("is the same as 'all' for copies with the same surroundings and orientation", async () => {
    // two copies 200 Å apart: a translation, no contact
    const builder = Structure.Builder({ label: 'two copies' });
    for (const shift of [0, 200]) {
      const op = SymmetryOperator.create(`t${shift}`, Mat4.fromTranslation(Mat4(), Vec3.create(shift, 0, 0)));
      for (const unit of model.units) builder.addWithOperator(unit, op);
    }
    const two = builder.getStructure();
    expect(two.unitSymmetryGroups.length).toBe(model.units.length);
    const all = await AccessibleSurfaceArea.compute(two, { symmetry: 'all' }).run();
    const perGroup = await AccessibleSurfaceArea.compute(two, { symmetry: 'per-group' }).run();
    expect(summary(perGroup)).toEqual(summary(all));
  });

  /**
   * The sphere points do not rotate with a copy, so 'all' differs between the copies of a closed assembly by the
   * sampling noise (most residues, up to ~15 Å² at 92 points). Rotating the whole assembly so that copy k sits in the
   * place of the representative gives copy k, under 'all', exactly the representative's surroundings and orientation:
   * what 'per-group' gives every copy. This holds because the group is closed (the rotated set is the same set); the
   * open half below fails it.
   */
  it("gives each copy of a closed point-group assembly the areas 'all' gives it in the representative's place", async () => {
    const n = model.units.length;
    const copies = [0, 1, 2, 3];
    const c4 = cyclicAssembly(model, 4, 20, copies);
    expect(c4.unitSymmetryGroups.length).toBe(n);
    const all = await AccessibleSurfaceArea.compute(c4, { symmetry: 'all' }).run();
    const perGroup = areasByUnit(c4, await AccessibleSurfaceArea.compute(c4, { symmetry: 'per-group' }).run());

    // not vacuous: the copies bury part of each other
    const single = await AccessibleSurfaceArea.compute(model).run();
    expect(summary(all).total).toBeLessThan(4 * summary(single).total - 500);

    // the copy of the protein (unit 0 of the model) that 'per-group' computed
    const r = representativeCopy(c4, areasByUnit(c4, all), perGroup, n, copies.length);
    for (let k = 0; k < 4; ++k) {
      const rotated = cyclicAssembly(
        model,
        4,
        20,
        copies.map((j) => j + r - k),
      );
      const inPlace = areasByUnit(rotated, await AccessibleSurfaceArea.compute(rotated, { symmetry: 'all' }).run());
      for (let u = k * n; u < (k + 1) * n; ++u) expect(perGroup.get(u)).toEqual(inPlace.get(u));
    }

    // an open assembly (two copies of the four) fails the same check: the copies' surroundings differ
    const open = cyclicAssembly(model, 4, 20, [0, 1]);
    const openAll = areasByUnit(open, await AccessibleSurfaceArea.compute(open, { symmetry: 'all' }).run());
    const openGroup = areasByUnit(open, await AccessibleSurfaceArea.compute(open, { symmetry: 'per-group' }).run());
    const or = representativeCopy(open, openAll, openGroup, n, 2);
    const k = 1 - or;
    const openRotated = cyclicAssembly(
      model,
      4,
      20,
      [0, 1].map((j) => j + or - k),
    );
    const openInPlace = areasByUnit(
      openRotated,
      await AccessibleSurfaceArea.compute(openRotated, { symmetry: 'all' }).run(),
    );
    expect(openGroup.get(k * n)).not.toEqual(openInPlace.get(k * n));
  });

  it('is exact for the most central unit of each group of a crystal block and copies it to the others', async () => {
    const all = areasByUnit(cells, await AccessibleSurfaceArea.compute(cells, { symmetry: 'all' }).run());
    const perGroup = areasByUnit(cells, await AccessibleSurfaceArea.compute(cells, { symmetry: 'per-group' }).run());
    const center = cells.boundary.sphere.center;
    const distance = (unit: Unit) =>
      Vec3.distance(Vec3.transformMat4(Vec3(), unit.boundary.sphere.center, unit.conformation.operator.matrix), center);
    let differs = false;
    for (const { units } of cells.unitSymmetryGroups) {
      const central = units.reduce((a, b) => (distance(b) < distance(a) ? b : a));
      const values = perGroup.get(central.id)!;
      expect(values).toEqual(all.get(central.id));
      for (const unit of units) {
        expect(perGroup.get(unit.id)).toEqual(values);
        if (maxDifference(all.get(unit.id)!, values) > 0) differs = true;
      }
    }
    // the approximation: copies at the edge of the block have other neighbours than the central one
    expect(differs).toBe(true);
  });

  it('recomputes when the parameter changes', async () => {
    const ctx = { runtime: SyncRuntimeContext, assetManager: new AssetManager() };
    const structure = await StructureSymmetry.buildSymmetryRange(
      model,
      Vec3.create(0, 0, 0),
      Vec3.create(1, 0, 0),
    ).run();
    await AccessibleSurfaceAreaProvider.attach(ctx, structure, { symmetry: 'per-group' });
    const perGroup = AccessibleSurfaceAreaProvider.get(structure).value!;
    await AccessibleSurfaceAreaProvider.attach(ctx, structure, { symmetry: 'all' });
    const all = AccessibleSurfaceAreaProvider.get(structure).value!;
    expect(all).not.toBe(perGroup);
    expect(summary(all)).toEqual({ residues: 1904, total: 56345, hash: 'f2c2e934' });
    expect(summary(perGroup)).not.toEqual({ residues: 1904, total: 56345, hash: 'f2c2e934' });
  });
});
