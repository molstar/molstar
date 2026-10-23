/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { parseLammpsTrajectory } from '@molstar/io/reader/lammps/traj/parser';
import { trajectoryFromLammpsTrajectory } from '@molstar/model/formats/structure/lammps-trajectory';
import type { Trajectory } from '@molstar/model/model/structure';
import { Model } from '@molstar/model/model/structure/model';
import { ModelPeriodicity } from '@molstar/model/model/structure/model/periodic';

// LAMMPS image flags (model.md §4 "Periodic frames") over data/examples/trajectory/periodic/
// long-chain.lammpstrj (its generator next to it): a 12-bead chain whose atom ids interleave it, so atoms with
// consecutive ids are more than half the box apart, drifting 6 Å a frame across the box's faces.

const text = fs.readFileSync(
  path.resolve(import.meta.dirname, '../../../../../../data/examples/trajectory/periodic/long-chain.lammpstrj'),
  'utf8',
);
/** The same dump without its image flags. */
const withoutImages = text.replace(/ ix iy iz/g, '').replace(/^(\d+ 1 1 \S+ \S+ \S+) -?\d+ -?\d+ -?\d+$/gm, '$1');

/** The chain's beads in id order (the atoms' order in the model): the generator's. */
const IdOrder = [0, 7, 1, 8, 2, 9, 3, 10, 4, 11, 5, 6];
const atomOfBead = (b: number) => IdOrder.indexOf(b);

async function trajectory(dump: string): Promise<Trajectory> {
  const parsed = await parseLammpsTrajectory(dump).run();
  if (parsed.isError) throw new Error(parsed.message);
  return trajectoryFromLammpsTrajectory(parsed.result).run();
}

/** The longest distance between consecutive beads of the chain. */
function longestBond(m: Model) {
  const { x, y, z } = m.atomicConformation;
  let max = 0;
  for (let b = 1; b < IdOrder.length; b++) {
    const i = atomOfBead(b - 1),
      j = atomOfBead(b);
    max = Math.max(max, Math.hypot(x[i] - x[j], y[i] - y[j], z[i] - z[j]));
  }
  return max;
}

const meanX = (m: Model) => Array.from(m.atomicConformation.x).reduce((a, v) => a + v, 0) / IdOrder.length;

describe('LAMMPS image flags', () => {
  it('make the chain whole in every frame, where the residue order cannot', async () => {
    const flagged = await trajectory(text);
    const plain = await trajectory(withoutImages);
    let broken = 0;
    for (let i = 0; i < flagged.frameCount; i++) {
      const frame = flagged.getFrameAtIndex(i) as Model;
      expect(Model.UnwrappedConformation.get(frame)).toBeDefined();
      // the bonds are 1.6 Å along x, 1 Å across: 1.89 Å
      expect(longestBond(ModelPeriodicity.makeWhole(frame))).toBeCloseTo(Math.hypot(1.6, 1), 3);
      const unflagged = plain.getFrameAtIndex(i) as Model;
      expect(Model.UnwrappedConformation.get(unflagged)).toBeUndefined();
      if (longestBond(ModelPeriodicity.makeWhole(unflagged)) > 5) broken++;
    }
    // the fixture needs them
    expect(broken).toBeGreaterThan(0);
  });

  it('leave the file as it is when not asked: the frame as read is wrapped', async () => {
    const frame = (await trajectory(text)).getFrameAtIndex(0) as Model;
    const { x } = frame.atomicConformation;
    // bead 0 at x = -12, read wrapped at 8
    expect(x[atomOfBead(0)]).toBeCloseTo(8, 4);
    expect(longestBond(frame)).toBeGreaterThan(10);
  });

  it('keep the chain from jumping, any frame first', async () => {
    const t = await trajectory(text);
    const inOrder: number[] = [];
    for (let i = 0; i < t.frameCount; i++) inOrder.push(meanX(await ModelPeriodicity.applyNoJump(t, i)));
    for (let i = 1; i < t.frameCount; i++) expect(inOrder[i] - inOrder[i - 1]).toBeCloseTo(6, 3);
    // a trajectory read again, frame 3 first: the same, with nothing walked
    const again = await trajectory(text);
    expect(meanX(await ModelPeriodicity.applyNoJump(again, 3))).toBeCloseTo(inOrder[3], 3);
    expect(longestBond(await ModelPeriodicity.applyNoJump(again, 3))).toBeLessThan(2);
  });
});
