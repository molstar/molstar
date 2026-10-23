/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { PeriodicBox } from '@molstar/core/math/geometry/periodic-box';
import { Cell } from '@molstar/core/math/geometry/spacegroup/cell';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { parseGRO } from '@molstar/io/reader/gro/parser';
import { trajectoryFromGRO } from '@molstar/model/formats/structure/gro';
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';
import { Coordinates, type Frame, Time } from '@molstar/model/model/structure/coordinates';
import { Model } from '../model.js';
import { ModelPeriodicity } from '../periodic.js';

// The fixture: data/examples/trajectory/periodic/split-peptide.gro (its generator next to it), a four-residue peptide
// and three waters in a rhombic dodecahedron, every atom wrapped into the box: split in most frames, the first among them.
const file = path.resolve(
  import.meta.dirname,
  '../../../../../../../data/examples/trajectory/periodic/split-peptide.gro',
);
const parsed = await parseGRO(fs.readFileSync(file, 'utf8')).run();
if (parsed.isError) throw new Error(parsed.message);
const trajectory = await trajectoryFromGRO(parsed.result).run();
const frames = Array.from({ length: trajectory.frameCount }, (_, i) => trajectory.getFrameAtIndex(i) as Model);

function distance(m: Model, i: number, j: number) {
  const { x, y, z } = m.atomicConformation;
  return Math.hypot(x[i] - x[j], y[i] - y[j], z[i] - z[j]);
}

/** The longest hop between an atom and its parent in the hierarchy plan: bonds and consecutive CA. */
function longestHop(m: Model) {
  const { parent } = ModelPeriodicity.planFromHierarchy(m.atomicHierarchy);
  let max = 0;
  for (let i = 0; i < parent.length; i++) if (parent[i] >= 0) max = Math.max(max, distance(m, i, parent[i]));
  return max;
}

describe('ModelPeriodicity', () => {
  it('reads the triclinic box of a GRO frame', () => {
    const box = Model.PeriodicBox.get(frames[0])!;
    expect(box).toBeDefined();
    // v3 = (15, 15, 21.2) Å: the dodecahedron, not a rectangular box
    expect(box.basis[6]).toBeCloseTo(15);
    expect(box.basis[7]).toBeCloseTo(15);
    expect(box.basis[8]).toBeCloseTo(21.2132, 3);
  });

  it('makes every frame whole', () => {
    let split = 0;
    for (const frame of frames) {
      if (longestHop(frame) > 5) split++;
      const whole = ModelPeriodicity.makeWhole(frame);
      expect(whole).not.toBe(frame);
      // the hops of the plan are bonds (< 1.5 Å) and peptide N–N links (3.75 Å)
      expect(longestHop(whole)).toBeLessThan(4);
    }
    // the fixture shows the problem
    expect(split).toBeGreaterThan(frames.length / 2);
    expect(longestHop(frames[0])).toBeGreaterThan(5);
  });

  it('is the same model for the same frame, and idempotent', () => {
    const whole = ModelPeriodicity.makeWhole(frames[0]);
    expect(ModelPeriodicity.makeWhole(frames[0])).toBe(whole);
    const again = ModelPeriodicity.applyWhole(
      ModelPeriodicity.getPlan(whole),
      Model.PeriodicBox.get(whole)!,
      whole.atomicConformation,
    );
    expect(Array.from(again.x)).toEqual(Array.from(whole.atomicConformation.x));
  });

  it('does not depend on which image an atom other than a root is in', () => {
    const frame = frames[0];
    const plan = ModelPeriodicity.getPlan(frame);
    const box = Model.PeriodicBox.get(frame)!;
    const { x, y, z } = frame.atomicConformation;
    // atom 7 (residue 2's CA) moved by a + c
    const moved = {
      ...frame.atomicConformation,
      x: Float32Array.from(x),
      y: Float32Array.from(y),
      z: Float32Array.from(z),
    };
    moved.x[7] += box.basis[0] + box.basis[6];
    moved.y[7] += box.basis[1] + box.basis[7];
    moved.z[7] += box.basis[2] + box.basis[8];
    const a = ModelPeriodicity.applyWhole(plan, box, frame.atomicConformation);
    const b = ModelPeriodicity.applyWhole(plan, box, moved);
    for (let i = 0; i < a.x.length; i++) {
      expect(b.x[i]).toBeCloseTo(a.x[i], 3);
      expect(b.y[i]).toBeCloseTo(a.y[i], 3);
      expect(b.z[i]).toBeCloseTo(a.z[i], 3);
    }
  });

  it('finds the polymer gaps again on whole coordinates', () => {
    // the topology frame is split: its ranges have a gap where the peptide crosses the face
    expect(frames[0].atomicRanges.gapRanges.length).toBeGreaterThan(0);
    const whole = ModelPeriodicity.makeWhole(frames[0]);
    expect(whole.atomicRanges.gapRanges.length).toBe(0);
    expect(whole.atomicRanges.polymerRanges.length).toBe(2);
  });

  it('agrees with a plan from topology bonds', () => {
    const frame = frames[0];
    const hierarchyPlan = ModelPeriodicity.planFromHierarchy(frame.atomicHierarchy);
    // bonds along the hierarchy plan's links, in reverse order: the BFS finds its own roots and order
    const a: number[] = [];
    const b: number[] = [];
    for (let i = hierarchyPlan.parent.length - 1; i >= 0; i--) {
      if (hierarchyPlan.parent[i] < 0) continue;
      a.push(i);
      b.push(hierarchyPlan.parent[i]);
    }
    const bonds = IndexPairBonds.fromData({
      pairs: { indexA: Column.ofIntArray(a), indexB: Column.ofIntArray(b) },
      count: frame.atomicHierarchy.atoms._rowCount,
    });
    const box = Model.PeriodicBox.get(frame)!;
    const fromBonds = ModelPeriodicity.applyWhole(
      ModelPeriodicity.planFromBonds(bonds, frame.atomicHierarchy.atoms._rowCount),
      box,
      frame.atomicConformation,
    );
    const fromHierarchy = ModelPeriodicity.applyWhole(hierarchyPlan, box, frame.atomicConformation);
    for (let i = 0; i < fromBonds.x.length; i++) {
      expect(fromBonds.x[i]).toBeCloseTo(fromHierarchy.x[i], 3);
      expect(fromBonds.y[i]).toBeCloseTo(fromHierarchy.y[i], 3);
      expect(fromBonds.z[i]).toBeCloseTo(fromHierarchy.z[i], 3);
    }
  });

  it('gives the frames of coordinates with a cell a box, and leaves the others alone', () => {
    const topology = frames[0];
    const n = topology.atomicHierarchy.atoms._rowCount;
    const { x, y, z } = topology.atomicConformation;
    const cell = Cell.fromBasis(Vec3.create(30, 0, 0), Vec3.create(0, 30, 0), Vec3.create(15, 15, 21.2132));
    const frame = (withCell: boolean): Frame => ({
      elementCount: n,
      time: Time(0, 'step'),
      x: Float32Array.from(x),
      y: Float32Array.from(y),
      z: Float32Array.from(z),
      ...(withCell ? { cell } : {}),
      xyzOrdering: { isIdentity: true },
    });
    const withBox = Model.trajectoryFromModelAndCoordinates(
      topology,
      Coordinates.create([frame(true)], Time(1, 'step'), Time(0, 'step')),
    ).getFrameAtIndex(0) as Model;
    expect(ModelPeriodicity.isApplicable(withBox)).toBe(true);
    expect(longestHop(ModelPeriodicity.makeWhole(withBox))).toBeLessThan(4);

    const without = Model.trajectoryFromModelAndCoordinates(
      topology,
      Coordinates.create([frame(false)], Time(1, 'step'), Time(0, 'step')),
    ).getFrameAtIndex(0) as Model;
    // the topology's symmetry is kept, but it is not a simulation box
    expect(ModelPeriodicity.isApplicable(without)).toBe(false);
    expect(ModelPeriodicity.makeWhole(without)).toBe(without);
  });

  it('leaves a frame in a placeholder box alone', () => {
    const frame = Model.withAtomicConformation(frames[0], frames[0].atomicConformation);
    Model.PeriodicBox.set(frame, PeriodicBox.create(Vec3.create(1, 0, 0), Vec3.create(0, 1, 0), Vec3.create(0, 0, 1)));
    expect(ModelPeriodicity.isApplicable(frame)).toBe(false);
    expect(ModelPeriodicity.makeWhole(frame)).toBe(frame);
  });
});
