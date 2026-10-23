/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, type Unit } from '@molstar/model/model/structure';
import { StructureSymmetry } from '../symmetry.js';

describe('structure boundary', () => {
  it('is that of the moved unit spheres, without moving copies of them (perf-audit-2 05/P15)', async () => {
    const file = path.resolve(import.meta.dirname, '../../../../../../../data/examples/1cbs_updated.cif');
    const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
    const s = await StructureSymmetry.buildSymmetryRange(
      Structure.ofModel(model),
      Vec3.create(0, 0, 0),
      Vec3.create(1, 1, 0),
    ).run();
    expect(s.units.some((u) => !u.conformation.operator.isIdentity)).toBe(true);
    for (const u of s.units) void u.boundary; // the unit spheres, before counting

    const transform = vi.spyOn(Sphere3D, 'transform');
    const { box, sphere } = s.boundary;
    expect(transform).not.toHaveBeenCalled();
    transform.mockRestore();

    // the same as the helper given each moved unit sphere
    const reference = new BoundaryHelper('98');
    const moved = (u: Unit) => Sphere3D.transform(Sphere3D(), u.boundary.sphere, u.conformation.operator.matrix);
    reference.reset();
    for (const u of s.units) reference.includeSphere(moved(u));
    reference.finishedIncludeStep();
    for (const u of s.units) reference.radiusSphere(moved(u));
    const expected = reference.getSphere();
    expect(Vec3.exactEquals(sphere.center, expected.center)).toBe(true);
    expect(sphere.radius).toBe(expected.radius);
    expect(Vec3.exactEquals(box.min, reference.getBox().min)).toBe(true);
    expect(Vec3.exactEquals(box.max, reference.getBox().max)).toBe(true);
  });
});
