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
    // per unit a copy for each transform, of its invariant (instances of one group), placed by it; at the identity the
    // unit itself, its id kept
    expect(new Set(s.units.map((u) => u.invariantId)).size).toBe(s.units.length);
    for (const u of s.units) {
      const p = u.conformation.position(u.elements[0], Vec3());
      const copies = instances.units.filter((c) => c.invariantId === u.invariantId);
      expect(copies.length).toBe(transforms.length);
      for (const transform of transforms) {
        const at = Vec3.transformMat4(Vec3(), p, transform);
        const copy = copies.find((c) => Vec3.distance(c.conformation.position(c.elements[0], Vec3()), at) < 1e-4);
        expect(copy).toBeDefined();
        expect(instances.unitMap.get(copy!.id)).toBe(copy);
        if (Mat4.isIdentity(transform)) expect(copy).toBe(u);
      }
    }
    // without the identity: numbered from 0
    const shifted = Structure.instances(s, [shift(40), shift(80)]).units.map((u) => u.id);
    expect(shifted).toEqual(shifted.map((_, i) => i));
    expect(Structure.instances(s, [Mat4.identity()]).units).toEqual(s.units);
  });

  it("keeps the ids of a component's units at the identity, by which its loci are remapped to the root", async () => {
    const s = await structure();
    // a component of the second unit alone (as of its ligand): its id not 0
    const component = Structure.create([s.units[1]], { parent: s });
    const [unit] = component.units;
    expect(unit.id).not.toBe(0);
    const instances = Structure.instances(component, [
      Mat4.identity(),
      Mat4.fromTranslation(Mat4(), Vec3.create(40, 0, 0)),
    ]);
    expect(instances.units[0]).toBe(unit);
    expect(instances.units[1].id).toBe(unit.id + 1);
    // the root's unit of that id is the identity copy's
    expect(s.unitMap.get(instances.units[0].id)).toBe(s.units[1]);
  });
});
