/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { parseMol2 } from '@molstar/io/reader/mol2/parser';
import { trajectoryFromMol2 } from '@molstar/model/formats/structure/mol2';
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';
import { Coordinates, type Frame, Time } from '@molstar/model/model/structure/coordinates';
import { Model } from '../model.js';

describe('a coordinate trajectory of a model', () => {
  it("keeps the model's explicit bonds on every frame (perf-audit-2 04/B7)", async () => {
    const file = path.resolve(import.meta.dirname, '../../../../../../../data/examples/adrenalin.mol2');
    const parsed = await parseMol2(fs.readFileSync(file, 'utf8'), 'adrenalin').run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMol2(parsed.result).run()).representative;
    const bonds = IndexPairBonds.Provider.get(model);
    expect(bonds?.bonds.edgeCount).toBeGreaterThan(0);

    const n = model.atomicHierarchy.atoms._rowCount;
    const { x, y, z } = model.atomicConformation;
    const frames: Frame[] = [0, 1, 2].map((i) => ({
      elementCount: n,
      time: Time(i, 'step'),
      x: Float32Array.from(x, (v) => v + i),
      y,
      z,
      xyzOrdering: { isIdentity: true },
    }));
    const trajectory = Model.trajectoryFromModelAndCoordinates(
      model,
      Coordinates.create(frames, Time(1, 'step'), Time(0, 'step')),
    );
    expect(trajectory.frameCount).toBe(3);
    for (let i = 0; i < 3; i++) {
      expect(IndexPairBonds.Provider.get(trajectory.getFrameAtIndex(i) as Model)).toBe(bonds);
    }
  });
});
