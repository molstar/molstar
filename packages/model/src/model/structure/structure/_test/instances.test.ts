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
import { Structure } from '../structure.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');

async function structure() {
  const parsed = await CIF.parse(fs.readFileSync(path.join(examples, '1cbs_updated.cif'), 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
  return Structure.ofModel(model);
}

describe('Structure.instances', () => {
  it('numbers every copy, the identity among the transforms', async () => {
    const s = await structure();
    expect(s.units.length).toBeGreaterThan(1);
    const shift = (x: number) => Mat4.fromTranslation(Mat4(), Vec3.create(x, 0, 0));
    // the identity in the middle: its copies' ids are those of other copies (formerly two units of one id)
    const transforms = [shift(40), Mat4.identity(), shift(80)];
    const instances = Structure.instances(s, transforms);

    const ids = instances.units.map((u) => u.id);
    expect(new Set(ids).size).toBe(s.units.length * transforms.length);
    expect(ids).toEqual(ids.map((_, i) => i));
    // per unit its copies in the transforms' order, of its invariant (instances of one group), placed by them
    for (let i = 0; i < s.units.length; ++i) {
      const u = s.units[i];
      const p = u.conformation.position(u.elements[0], Vec3());
      for (let t = 0; t < transforms.length; ++t) {
        const copy = instances.units[i * transforms.length + t];
        expect(copy.invariantId).toBe(u.invariantId);
        const q = copy.conformation.position(copy.elements[0], Vec3());
        expect(Vec3.distance(q, Vec3.transformMat4(Vec3(), p, transforms[t]))).toBeLessThan(1e-4);
        expect(instances.unitMap.get(copy.id)).toBe(copy);
      }
    }
    // a copy at the identity of a unit whose id is the next: the unit itself
    expect(Structure.instances(s, [Mat4.identity()]).units).toEqual(s.units);
  });
});
