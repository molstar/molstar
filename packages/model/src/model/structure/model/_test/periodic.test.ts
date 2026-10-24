/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { setFlagsFromString } from 'node:v8';
import { runInNewContext } from 'node:vm';
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

  it('does not keep the frames it made once the trajectory is dropped (perf-audit-3 05/P23)', async () => {
    // a trajectory of its own, made whole and in the box frame by frame, then dropped
    const made = await (async () => {
      const own = await trajectoryFromGRO(parsed.result).run();
      const refs: WeakRef<Model>[] = [];
      for (let i = 0; i < own.frameCount; i++) {
        const frame = own.getFrameAtIndex(i) as Model;
        refs.push(new WeakRef(ModelPeriodicity.apply(frame, ModelPeriodicity.InBox)));
        // still the same model while the topology is alive
        expect(ModelPeriodicity.apply(frame, ModelPeriodicity.InBox)).toBe(refs[i].deref());
      }
      return refs;
    })();
    setFlagsFromString('--expose-gc');
    const gc: () => void = runInNewContext('gc');
    // WeakRef targets are kept until the job that made them ends
    await new Promise((r) => setTimeout(r, 0));
    gc();
    await new Promise((r) => setTimeout(r, 0));
    gc();
    expect(made.filter((r) => r.deref() !== undefined).length).toBe(0);
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

  /** Each molecule's center (the mean of its atoms) in fractional coordinates of the frame's box. */
  function moleculeCenters(m: Model) {
    const plan = ModelPeriodicity.getPlan(m);
    const box = Model.PeriodicBox.get(m)!;
    const { x, y, z } = m.atomicConformation;
    const sums = Array.from({ length: plan.moleculeCount }, () => [0, 0, 0, 0]);
    for (let i = 0; i < x.length; i++) {
      const s = sums[plan.molecule[i]];
      s[0] += x[i];
      s[1] += y[i];
      s[2] += z[i];
      s[3]++;
    }
    return sums.map(([sx, sy, sz, n]) => Vec3.transformMat3(Vec3(), Vec3.create(sx / n, sy / n, sz / n), box.inverse));
  }

  it('puts every molecule in the box', () => {
    for (const frame of frames) {
      const inBox = ModelPeriodicity.apply(frame, ModelPeriodicity.InBox);
      expect(longestHop(inBox)).toBeLessThan(4);
      // the peptide and three waters
      const centers = moleculeCenters(inBox);
      expect(centers).toHaveLength(4);
      for (const f of centers) for (const v of f) expect(v >= 0 && v < 1).toBe(true);
    }
  });

  it('centers the box on a selection', () => {
    // the peptide's atoms (residues 1 to 4, five atoms each)
    const center: ModelPeriodicity.Center = { atoms: Array.from({ length: 20 }, (_, i) => i) };
    for (const frame of frames) {
      const centered = ModelPeriodicity.apply(frame, { kind: 'box', center });
      expect(ModelPeriodicity.apply(frame, { kind: 'box', center })).toBe(centered);
      expect(longestHop(centered)).toBeLessThan(4);
      const centers = moleculeCenters(centered);
      // the peptide's mean at the box's center, the waters in the box around it
      for (const v of centers[0]) expect(v).toBeCloseTo(0.5, 5);
      for (const f of centers.slice(1)) for (const v of f) expect(v >= 0 && v < 1).toBe(true);
    }
  });

  /** The peptide's mean (its first 20 atoms), Å. */
  function peptideCenter(m: Model) {
    const { x, y, z } = m.atomicConformation;
    const c = Vec3();
    for (let i = 0; i < 20; i++) Vec3.add(c, c, Vec3.create(x[i], y[i], z[i]));
    return Vec3.scale(c, c, 1 / 20);
  }

  it('keeps molecules from jumping between frames', async () => {
    const steps = (made: Model[]) =>
      made.slice(1).map((m, i) => Vec3.distance(peptideCenter(m), peptideCenter(made[i])));
    // made whole, the peptide jumps a box length where it crosses a face; without jumps it drifts 4.5 Å a frame
    const whole = frames.map((f) => ModelPeriodicity.makeWhole(f));
    expect(Math.max(...steps(whole))).toBeGreaterThan(15);
    const sequential: Model[] = [];
    for (let i = 0; i < frames.length; i++) sequential.push(await ModelPeriodicity.applyNoJump(trajectory, i));
    expect(Math.max(...steps(sequential))).toBeLessThan(8);
    for (const m of sequential) expect(longestHop(m)).toBeLessThan(4);
    // frame 0 is as made whole
    expect(Array.from(sequential[0].atomicConformation.x)).toEqual(Array.from(whole[0].atomicConformation.x));

    // a trajectory read again: frame 7 first, then back to frame 3, the same as walked in order
    const again = await trajectoryFromGRO(parsed.result).run();
    for (const i of [7, 3, 5]) {
      const m = await ModelPeriodicity.applyNoJump(again, i);
      expect(Vec3.distance(peptideCenter(m), peptideCenter(sequential[i]))).toBeLessThan(1e-3);
    }
  });

  it('seeks without jumps past its checkpoints', async () => {
    // 150 frames of the first one drifting 3 Å a frame along a, every atom wrapped into the box as an engine writes it
    const topology = frames[0];
    const box = Model.PeriodicBox.get(topology)!;
    const cell = Cell.fromBasis(
      Vec3.create(box.basis[0], box.basis[1], box.basis[2]),
      Vec3.create(box.basis[3], box.basis[4], box.basis[5]),
      Vec3.create(box.basis[6], box.basis[7], box.basis[8]),
    );
    const whole = ModelPeriodicity.makeWhole(topology).atomicConformation;
    const n = whole.x.length;
    const drifting = Array.from({ length: 150 }, (_, t): Frame => {
      const x = new Float32Array(n);
      const y = new Float32Array(n);
      const z = new Float32Array(n);
      const p = Vec3();
      for (let i = 0; i < n; i++) {
        Vec3.set(p, whole.x[i] + 3 * t, whole.y[i], whole.z[i]);
        Vec3.transformMat3(p, p, box.inverse);
        Vec3.set(p, p[0] - Math.floor(p[0]), p[1] - Math.floor(p[1]), p[2] - Math.floor(p[2]));
        Vec3.transformMat3(p, p, box.basis);
        [x[i], y[i], z[i]] = [p[0], p[1], p[2]];
      }
      return { elementCount: n, time: Time(t, 'step'), x, y, z, cell, xyzOrdering: { isIdentity: true } };
    });
    const coordinates = () => Coordinates.create(drifting, Time(1, 'step'), Time(0, 'step'));
    const walked = Model.trajectoryFromModelAndCoordinates(topology, coordinates());
    const inOrder: Vec3[] = [];
    for (let i = 0; i < 150; i++) inOrder.push(peptideCenter(await ModelPeriodicity.applyNoJump(walked, i)));
    // 3 Å a frame along x, no jumps
    for (let i = 1; i < 150; i++) expect(inOrder[i][0] - inOrder[i - 1][0]).toBeCloseTo(3, 3);

    const seeked = Model.trajectoryFromModelAndCoordinates(topology, coordinates());
    for (const i of [140, 70, 130, 65, 0, 149]) {
      const c = peptideCenter(await ModelPeriodicity.applyNoJump(seeked, i));
      expect(Vec3.distance(c, inOrder[i])).toBeLessThan(1e-2);
    }
  });

  describe('distance links across the box', () => {
    // data/examples/trajectory/periodic/split-ligand.gro (its generator next to it): a ligand of two residues whose
    // joining bond (LG1 C4, atom 4 - LG2 C1, atom 5) crosses a face in the first frame, a water (9-11) and a sodium (12)
    const ligandFile = path.resolve(
      import.meta.dirname,
      '../../../../../../../data/examples/trajectory/periodic/split-ligand.gro',
    );
    const ligandFrames = async () => {
      const p = await parseGRO(fs.readFileSync(ligandFile, 'utf8')).run();
      if (p.isError) throw new Error(p.message);
      const t = await trajectoryFromGRO(p.result).run();
      return Array.from({ length: t.frameCount }, (_, i) => t.getFrameAtIndex(i) as Model);
    };

    it('join the residues the hierarchy leaves apart, not the water or the ion', async () => {
      const [first] = await ligandFrames();
      const plan = ModelPeriodicity.getPlan(first);
      expect(plan.molecule[4]).toBe(plan.molecule[5]);
      expect(plan.molecule[9]).not.toBe(plan.molecule[0]);
      expect(plan.molecule[12]).not.toBe(plan.molecule[0]);
      expect(plan.molecule[12]).not.toBe(plan.molecule[9]);
      // the hierarchy alone: two molecules for the ligand
      const alone = ModelPeriodicity.planFromHierarchy(first.atomicHierarchy);
      expect(alone.molecule[4]).not.toBe(alone.molecule[5]);
    });

    it('make the ligand whole in every frame', async () => {
      const frames = await ligandFrames();
      const joint = (m: Model) => distance(m, 4, 5);
      expect(joint(frames[0])).toBeGreaterThan(10);
      for (const f of frames) expect(joint(ModelPeriodicity.makeWhole(f))).toBeCloseTo(1.487, 2);
      // the hierarchy alone leaves the first frame split
      const box = Model.PeriodicBox.get(frames[0])!;
      const alone = ModelPeriodicity.applyWhole(
        ModelPeriodicity.planFromHierarchy(frames[0].atomicHierarchy),
        box,
        frames[0].atomicConformation,
      );
      const { x, y, z } = alone;
      expect(Math.hypot(x[4] - x[5], y[4] - y[5], z[4] - z[5])).toBeGreaterThan(10);
    });

    it('leave the ion to move into the box on its own', async () => {
      const frames = await ligandFrames();
      for (const f of frames) {
        const inBox = ModelPeriodicity.apply(f, ModelPeriodicity.InBox);
        const box = Model.PeriodicBox.get(inBox)!;
        const { x, y, z } = inBox.atomicConformation;
        const p = Vec3.transformMat3(Vec3(), Vec3.create(x[12], y[12], z[12]), box.inverse);
        for (const v of p) expect(v >= 0 && v < 1).toBe(true);
        expect(distance(inBox, 4, 5)).toBeCloseTo(1.487, 2);
      }
    });
  });
});
