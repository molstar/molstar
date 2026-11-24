/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '../../structure.js';
import { StructureSymmetry } from '../../symmetry.js';
import { StructureUnitTransforms } from '../unit-transforms.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../../data/examples');

describe('StructureUnitTransforms', () => {
  it("gives the transforms of groups other than the structure's own, by their units (a representation's groups)", async () => {
    const parsed = await CIF.parse(fs.readFileSync(path.join(examples, '1cbs_updated.cif'), 'utf8')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
    const structure = Structure.ofModel(model);
    const transforms = new StructureUnitTransforms(structure);
    const unit = structure.units[1];
    const m = Mat4.fromTranslation(Mat4(), Vec3.create(1, 2, 3));
    transforms.setTransform(m, unit);

    // new group objects, as the units representations compute them
    const groups = StructureSymmetry.computeChildAwareTransformGroups(structure);
    expect(groups.some((g) => structure.unitSymmetryGroups.includes(g))).toBe(false);
    const group = groups.find((g) => g.units.includes(unit))!;
    const j = group.units.indexOf(unit);
    const array = transforms.getSymmetryGroupTransforms(group);
    expect(Array.from(array.subarray(j * 16, j * 16 + 16))).toEqual(Array.from(m));
    // the structure's own group: the same transforms
    const own = structure.unitSymmetryGroups.find((g) => g.units.includes(unit))!;
    expect(Array.from(transforms.getSymmetryGroupTransforms(own))).toEqual(Array.from(array));

    // changed later: read again
    Mat4.fromTranslation(m, Vec3.create(4, 5, 6));
    transforms.setTransform(m, unit);
    const again = transforms.getSymmetryGroupTransforms(group);
    expect(Array.from(again.subarray(j * 16, j * 16 + 16))).toEqual(Array.from(m));
  });
});
