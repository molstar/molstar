/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { parseLammpsData } from '@molstar/io/reader/lammps/data/parser';
import { parseLammpsTrajectory } from '@molstar/io/reader/lammps/traj/parser';
import { trajectoryFromLammpsData } from '@molstar/model/formats/structure/lammps-data';
import { trajectoryFromLammpsTrajectory } from '@molstar/model/formats/structure/lammps-trajectory';
import { Model } from '@molstar/model/model/structure/model';
import { ModelPeriodicity } from '@molstar/model/model/structure/model/periodic';

// A triclinic box from 10 to 30 in x, -5 to 15 in y, 0 to 20 in z, tilts xy = 4, xz = -2, yz = 3. The dump's bound
// lines are of the bounding box: xlo - min(0, xy, xz, xy + xz) = 10 - (-2), xhi - max(0, ...) = 30 + 4, ylo, yhi by yz.
function dump(flags: string, atoms: string) {
  return `ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS xy xz yz ${flags}
8.0 34.0 4.0
-5.0 18.0 -2.0
0.0 20.0 3.0
ITEM: ATOMS id type ${atoms}
1 1 0.5 0.5 0.5
2 1 0.25 0 0
`;
}

async function trajectory(text: string) {
  const parsed = await parseLammpsTrajectory(text).run();
  if (parsed.isError) throw new Error(parsed.message);
  return { file: parsed.result, trajectory: await trajectoryFromLammpsTrajectory(parsed.result).run() };
}

describe('LAMMPS boxes', () => {
  it('reads a triclinic dump box: its true bounds and tilt', async () => {
    const { file } = await trajectory(dump('pp pp pp', 'xs ys zs'));
    const [box] = file.bounds;
    expect(box.tilt).toEqual([4, -2, 3]);
    expect(box.lower).toEqual([10, -5, 0]);
    expect(box.length).toEqual([20, 20, 20]);
    expect(box.periodicity).toEqual(['pp', 'pp', 'pp']);
  });

  it('gives the frames the box at its corner, and scaled coordinates through the tilt', async () => {
    const { trajectory: t } = await trajectory(dump('pp pp pp', 'xs ys zs'));
    const frame = t.getFrameAtIndex(0) as Model;
    const box = Model.PeriodicBox.get(frame)!;
    expect(Array.from(box.origin)).toEqual([10, -5, 0]);
    // a = (20, 0, 0), b = (4, 20, 0), c = (-2, 3, 20), as columns
    expect(Array.from(box.basis)).toEqual([20, 0, 0, 4, 20, 0, -2, 3, 20]);
    // (0.5, 0.5, 0.5): the corner plus half of a + b + c
    const { x, y, z } = frame.atomicConformation;
    expect([x[0], y[0], z[0]]).toEqual([10 + 11, -5 + 11.5, 10]);
    expect([x[1], y[1], z[1]]).toEqual([15, -5, 0]);

    // in the box, by its corner: both atoms' fractional coordinates in [0, 1)
    const inBox = ModelPeriodicity.apply(frame, ModelPeriodicity.InBox);
    expect(inBox.atomicConformation.x[0]).toBeCloseTo(21);
  });

  it('has no box when a dimension is not periodic', async () => {
    const { trajectory: t } = await trajectory(dump('pp pp ff', 'x y z'));
    expect(Model.PeriodicBox.get(t.getFrameAtIndex(0) as Model)).toBeUndefined();
  });

  it('reads a data file box with tilt factors', async () => {
    const parsed = await parseLammpsData(`LAMMPS data file

1 atoms
1 atom types

10.0 30.0 xlo xhi
-5.0 15.0 ylo yhi
0.0 20.0 zlo zhi
4.0 -2.0 3.0 xy xz yz

Atoms # atomic

1 1 12.0 0.0 1.0
`).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromLammpsData(parsed.result).run()).representative;
    const box = Model.PeriodicBox.get(model)!;
    expect(Array.from(box.origin)).toEqual([10, -5, 0]);
    expect(Array.from(box.basis)).toEqual([20, 0, 0, 4, 20, 0, -2, 3, 20]);
  });
});
