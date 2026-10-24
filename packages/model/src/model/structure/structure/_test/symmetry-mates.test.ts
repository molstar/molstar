/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { SymmetryOperator } from '@molstar/core/math/geometry/symmetry-operator';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '../structure.js';
import { StructureSymmetry } from '../symmetry.js';
import type { Unit } from '../unit.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');

const key = (invariantId: number, m: Mat4) => `${invariantId}|${Array.from(m, (v) => v.toFixed(3)).join(',')}`;

/** The mates by their definition: a unit under an operator is kept if any atom lies within its sphere radius + `radius` */
function expectedMates(s: Structure, operators: Mat4[], radius: number) {
  const atoms: Vec3[] = [];
  for (const u of s.units) {
    for (let i = 0; i < u.elements.length; i++) atoms.push(u.conformation.position(u.elements[i], Vec3()));
  }
  const keys = new Set<string>();
  const c = Vec3(),
    m = Mat4();
  for (const op of operators) {
    for (const u of s.units) {
      const sphere = u.lookup3d.boundary.sphere;
      Mat4.mul(m, op, u.conformation.operator.matrix);
      Vec3.transformMat4(c, sphere.center, m);
      const r = sphere.radius + radius;
      if (atoms.some((a) => Vec3.distance(a, c) <= r)) keys.add(key(u.invariantId, m));
    }
  }
  return [...keys].sort();
}

const mateKeys = (s: Structure) =>
  [...new Set(s.units.map((u: Unit) => key(u.invariantId, u.conformation.operator.matrix)))].sort();

describe('symmetry mates', () => {
  it('a structure with non-identity unit operators gets the mates of its world positions (perf-audit-2 05/B8)', async () => {
    const parsed = await CIF.parse(fs.readFileSync(path.join(examples, '1cbs_updated.cif'), 'utf8')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
    const root = Structure.ofModel(model);

    // every operator the search tries: all of them pass at a huge radius
    const all = await StructureSymmetry.builderSymmetryMates(root, 1e6).run();
    const operators: Mat4[] = [];
    const seen = new Set<string>();
    for (const u of all.units) {
      const k = key(0, u.conformation.operator.matrix);
      if (!seen.has(k)) {
        seen.add(k);
        operators.push(u.conformation.operator.matrix);
      }
    }
    expect(operators.length).toBeGreaterThan(1);

    const radius = 5;
    const mates = await StructureSymmetry.builderSymmetryMates(root, radius).run();
    expect(mateKeys(mates)).toEqual(expectedMates(root, operators, radius));

    // the same structure moved by 60 Å along x: its units have a non-identity operator
    const t = SymmetryOperator.create('t', Mat4.fromTranslation(Mat4(), Vec3.create(60, 0, 0)));
    const moved = Structure.create(root.units.map((u) => u.applyOperator(u.id, t, true)));
    const movedMates = await StructureSymmetry.builderSymmetryMates(moved, radius).run();
    expect(mateKeys(movedMates)).toEqual(expectedMates(moved, operators, radius));
  });
});
