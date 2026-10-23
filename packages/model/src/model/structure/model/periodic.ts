/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PeriodicBox } from '@molstar/core/math/geometry/periodic-box';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { UUID } from '@molstar/core/util/uuid';
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';
import { type RuntimeContext, Task } from '@molstar/core/task';
import type { Trajectory } from '../trajectory.js';
import { Model } from './model.js';
import type { AtomicConformation, AtomicHierarchy, AtomicRanges } from './properties/atomic.js';
import { getAtomicRanges } from './properties/utils/atomic-ranges.js';
import { isPolymer } from './types.js';

// Molecules split by a simulation's periodic box made whole again (.v6/periodic-trajectories/plan.md §3): every atom
// moved by a lattice vector to the image nearest to its parent in a spanning forest of the topology, built once. Then,
// optionally, each whole molecule moved by a lattice vector into the box (plan.md §1, wrap), around a selection, or
// kept from jumping a box length between frames (no jump).

export { ModelPeriodicity };

/**
 * How to make a topology's molecules whole: each atom's parent, the atom it is imaged next to (-1: a root, it stays
 * where it is), and the order to visit them in, parents first.
 */
interface UnwrapPlan {
  readonly parent: Int32Array;
  /** Atoms in visiting order; absent: in index order (every parent comes before its children) */
  readonly order?: Int32Array;
  /** Each atom's molecule: the tree of the forest it is in */
  readonly molecule: Int32Array;
  readonly moleculeCount: number;
  /** The topology's ranges computed on a whole frame (polymer gaps): shared by its frames once made whole */
  ranges?: AtomicRanges;
}

namespace ModelPeriodicity {
  /** Atoms to center the box on (a selection's): its identity is the cache key, so one object per topology. */
  export interface Center {
    readonly atoms: ArrayLike<number>;
  }

  /**
   * `whole`: each molecule in one piece; `box`: also each molecule's center in the box, or, around a center, the
   * center's molecules kept together at the box's center and every other molecule's center in the box around it.
   */
  export type Mode = { readonly kind: 'whole' } | { readonly kind: 'box'; readonly center?: Center };
  export const Whole: Mode = { kind: 'whole' };
  export const InBox: Mode = { kind: 'box' };

  const plans = new WeakMap<object, UnwrapPlan>();
  // The frames made, by source frame and mode: the same model for the same frame while it is recent. Not every frame:
  // a trajectory read in full would be held twice.
  const CacheSize = 16;
  const recent = new Map<Model, Map<unknown, Model>>();

  function cached(model: Model, key: unknown) {
    const byMode = recent.get(model);
    if (!byMode) return;
    recent.delete(model);
    recent.set(model, byMode);
    return byMode.get(key);
  }

  function remember(model: Model, key: unknown, made: Model) {
    let byMode = recent.get(model);
    if (!byMode) {
      byMode = new Map();
      recent.set(model, byMode);
      if (recent.size > CacheSize) recent.delete(recent.keys().next().value!);
    }
    byMode.set(key, made);
  }

  /**
   * The shortest box vector a frame is made whole in, Å: the plan's hops (bonds, 3.8 Å between consecutive residues)
   * must be under half of it. A smaller box is a placeholder (a 1 Å cell), not a simulation's.
   */
  export const MinBoxVector = 10;

  /** The frame's box, when it can be made whole in it. */
  function boxOf(model: Model): PeriodicBox | undefined {
    const box = Model.PeriodicBox.get(model);
    if (!box || model.atomicHierarchy.atoms._rowCount === 0) return;
    const m = box.basis;
    const shortest = Math.min(Math.hypot(m[0], m[1], m[2]), Math.hypot(m[3], m[4], m[5]), Math.hypot(m[6], m[7], m[8]));
    return shortest >= MinBoxVector ? box : undefined;
  }

  /** A frame can be made whole: it has atoms and a simulation box of at least `MinBoxVector`. */
  export function isApplicable(model: Model): boolean {
    return !!boxOf(model);
  }

  /** The plan for the model's topology: by its topology bonds when it has them, else by its hierarchy. */
  export function getPlan(model: Model): UnwrapPlan {
    const bonds = IndexPairBonds.Provider.get(model);
    // frames of one topology share the hierarchy (and the bonds); the frames of a multi-model file have their own
    const key = bonds ?? model.atomicHierarchy;
    let plan = plans.get(key);
    if (!plan) {
      plan = bonds
        ? planFromBonds(bonds, model.atomicHierarchy.atoms._rowCount)
        : planFromHierarchy(model.atomicHierarchy);
      plans.set(key, plan);
    }
    return plan;
  }

  /**
   * The frame with its molecules whole: a new model (the same one for the same frame while it is recent), or the frame
   * itself when it has no simulation box (`isApplicable`).
   */
  export function makeWhole(model: Model): Model {
    return apply(model, Whole);
  }

  /** The frame made whole, and with `box` its molecules in the box; the frame itself when it has no simulation box. */
  export function apply(model: Model, mode: Mode): Model {
    const box = boxOf(model);
    if (!box) return model;
    const key = mode.kind === 'whole' ? 'whole' : (mode.center ?? 'box');
    const hit = cached(model, key);
    if (hit) return hit;

    const plan = getPlan(model);
    let made: Model;
    if (mode.kind === 'whole') {
      const conformation = applyWhole(plan, box, model.atomicConformation);
      // polymer gaps were found on the topology's coordinates, which may have been split: found again on whole ones
      if (!plan.ranges) {
        plan.ranges = getAtomicRanges(model.atomicHierarchy, model.entities, conformation, model.sequence);
      }
      made = Model.withAtomicConformation(model, conformation, plan.ranges);
    } else {
      const whole = apply(model, Whole);
      // whole molecules moved whole: the whole frame's ranges hold
      const conformation = applyBox(plan, box, whole.atomicConformation, mode.center?.atoms);
      made = Model.withAtomicConformation(whole, conformation);
    }
    remember(model, key, made);
    return made;
  }

  /**
   * Each molecule of whole coordinates moved by a lattice vector so its center (the mean of its atoms) is in the box,
   * in new arrays. With `center` atoms: their molecules kept together around the first one's, everything translated
   * so the atoms' mean is at the box's center, and every other molecule's center put in the box around it.
   */
  export function applyBox(
    plan: UnwrapPlan,
    box: PeriodicBox,
    whole: AtomicConformation,
    center?: ArrayLike<number>,
  ): AtomicConformation {
    const { molecule, moleculeCount } = plan;
    const x = Float32Array.from(whole.x);
    const y = Float32Array.from(whole.y);
    const z = Float32Array.from(whole.z);
    const n = x.length;

    // molecule centers, then each molecule's shift into the box
    const centers = moleculeCenters(plan, whole);
    const shift = new Float64Array(moleculeCount * 3);
    const { basis, inverse, origin } = box;
    const boxCenter = Vec3.add(Vec3(), origin, Vec3.transformMat3(Vec3(), Vec3.create(0.5, 0.5, 0.5), basis));
    const t = Vec3();

    const pinned = center && center.length ? new Uint8Array(moleculeCount) : undefined;
    if (pinned) {
      // the center's molecules at the images nearest to its first molecule: one piece (a complex split by the box)
      const first = molecule[center![0]];
      const ref = Vec3.fromArray(Vec3(), centers, 3 * first);
      const d = Vec3();
      for (let k = 0; k < center!.length; k++) {
        const m = molecule[center![k]];
        if (pinned[m]) continue;
        pinned[m] = 1;
        Vec3.sub(d, Vec3.fromArray(d, centers, 3 * m), ref);
        // the shift from the molecule's center to its image nearest the reference: d's minimum image minus d
        const raw = Vec3.clone(d);
        PeriodicBox.minimumImage(d, box, d);
        Vec3.toArray(Vec3.sub(d, d, raw), shift, 3 * m);
      }
      const sum = Vec3();
      for (let k = 0; k < center!.length; k++) {
        const i = center![k];
        const m = 3 * molecule[i];
        sum[0] += x[i] + shift[m];
        sum[1] += y[i] + shift[m + 1];
        sum[2] += z[i] + shift[m + 2];
      }
      Vec3.sub(t, boxCenter, Vec3.scale(sum, sum, 1 / center!.length));
    }

    const p = Vec3();
    for (let m = 0; m < moleculeCount; m++) {
      if (pinned?.[m]) continue;
      // fractional coordinates of the moved center, relative to the box's corner; whole box vectors back
      p[0] = centers[3 * m] + t[0] - origin[0];
      p[1] = centers[3 * m + 1] + t[1] - origin[1];
      p[2] = centers[3 * m + 2] + t[2] - origin[2];
      Vec3.transformMat3(p, p, inverse);
      p[0] = -Math.floor(p[0]);
      p[1] = -Math.floor(p[1]);
      p[2] = -Math.floor(p[2]);
      Vec3.transformMat3(p, p, basis);
      Vec3.toArray(p, shift, 3 * m);
    }

    for (let i = 0; i < n; i++) {
      const m = 3 * molecule[i];
      x[i] += t[0] + shift[m];
      y[i] += t[1] + shift[m + 1];
      z[i] += t[2] + shift[m + 2];
    }
    return { ...whole, id: UUID.create22(), x, y, z };
  }

  /** Each molecule's center, the mean of its atoms: x, y, z per molecule. */
  export function moleculeCenters(plan: UnwrapPlan, conformation: AtomicConformation): Float64Array {
    const { molecule, moleculeCount } = plan;
    const { x, y, z } = conformation;
    const centers = new Float64Array(moleculeCount * 3);
    const counts = new Int32Array(moleculeCount);
    for (let i = 0, n = molecule.length; i < n; i++) {
      const m = molecule[i] * 3;
      centers[m] += x[i];
      centers[m + 1] += y[i];
      centers[m + 2] += z[i];
      counts[molecule[i]]++;
    }
    for (let m = 0; m < moleculeCount; m++) {
      const c = counts[m] || 1;
      centers[3 * m] /= c;
      centers[3 * m + 1] /= c;
      centers[3 * m + 2] /= c;
    }
    return centers;
  }

  // No jump (plan.md §1): each molecule's lattice shift, in box vectors, chosen so its center moves less than half a
  // box from one frame to the next. Sequential by nature: frame k's shifts need frame k - 1's. Kept per trajectory:
  // the last frame's (a playing trajectory steps one frame, either way) and a checkpoint every `NoJumpEvery` frames
  // (a seek walks at most that far), not every frame's (a molecule count of shifts per frame).
  export const NoJumpEvery = 64;

  interface NoJumpState {
    readonly checkpoints: Map<number, Int32Array>;
    last?: {
      readonly index: number;
      readonly shifts: Int32Array;
      readonly centers: Float64Array;
      readonly box?: PeriodicBox;
    };
  }
  const noJumpStates = new WeakMap<Trajectory, NoJumpState>();

  /**
   * The trajectory's frame made whole, each molecule moved by whole box vectors so its center never jumps a box
   * length from the frame before: molecules leave the box as they diffuse. Frame 0 is as made whole. A frame without
   * a simulation box is returned as it is.
   */
  export async function applyNoJump(trajectory: Trajectory, index: number, ctx?: RuntimeContext): Promise<Model> {
    const frame = await frameAt(trajectory, index, ctx);
    const box = boxOf(frame);
    if (!box) return frame;
    const hit = cached(frame, 'nojump');
    if (hit) return hit;

    const shifts = await noJumpShifts(trajectory, index, ctx);
    const whole = apply(frame, Whole);
    const plan = getPlan(frame);
    const x = Float32Array.from(whole.atomicConformation.x);
    const y = Float32Array.from(whole.atomicConformation.y);
    const z = Float32Array.from(whole.atomicConformation.z);
    const b = box.basis;
    const { molecule } = plan;
    for (let i = 0, n = x.length; i < n; i++) {
      const m = 3 * molecule[i];
      const na = shifts[m],
        nb = shifts[m + 1],
        nc = shifts[m + 2];
      if (na === 0 && nb === 0 && nc === 0) continue;
      x[i] += na * b[0] + nb * b[3] + nc * b[6];
      y[i] += na * b[1] + nb * b[4] + nc * b[7];
      z[i] += na * b[2] + nb * b[5] + nc * b[8];
    }
    const made = Model.withAtomicConformation(whole, { ...whole.atomicConformation, id: UUID.create22(), x, y, z });
    remember(frame, 'nojump', made);
    return made;
  }

  async function frameAt(trajectory: Trajectory, index: number, ctx?: RuntimeContext) {
    const f = trajectory.getFrameAtIndex(index);
    return ctx ? Task.resolveInContext(f, ctx) : Task.is(f) ? f.run() : f;
  }

  /** The centers of the frame's whole molecules, and its box. */
  async function wholeCenters(trajectory: Trajectory, index: number, ctx?: RuntimeContext) {
    const frame = await frameAt(trajectory, index, ctx);
    const box = boxOf(frame);
    const whole = apply(frame, Whole);
    return { centers: moleculeCenters(getPlan(frame), whole.atomicConformation), box };
  }

  /** Frame `index`'s shifts: walked from the nearest frame known, the last one or a checkpoint. */
  async function noJumpShifts(trajectory: Trajectory, index: number, ctx?: RuntimeContext): Promise<Int32Array> {
    let state = noJumpStates.get(trajectory);
    if (!state) {
      state = { checkpoints: new Map() };
      noJumpStates.set(trajectory, state);
    }
    const { checkpoints } = state;
    if (state.last?.index === index) return state.last.shifts;
    const known = checkpoints.get(index);
    if (known) return known;

    // the nearest start: the last frame, a checkpoint, or frame 0 (no shifts)
    let start = 0;
    let shifts: Int32Array | undefined;
    let centers: Float64Array | undefined;
    let box: PeriodicBox | undefined;
    for (const c of checkpoints.keys()) if (Math.abs(c - index) < Math.abs(start - index)) start = c;
    const last = state.last;
    if (last && Math.abs(last.index - index) <= Math.abs(start - index)) {
      start = last.index;
      shifts = Int32Array.from(last.shifts);
      centers = last.centers;
      box = last.box;
    } else {
      const atStart = await wholeCenters(trajectory, start, ctx);
      centers = atStart.centers;
      box = atStart.box;
      shifts = start === 0 ? new Int32Array(centers.length) : Int32Array.from(checkpoints.get(start)!);
      if (start === 0) checkpoints.set(0, Int32Array.from(shifts));
    }

    const step = index > start ? 1 : -1;
    for (let k = start + step; step > 0 ? k <= index : k >= index; k += step) {
      const next = await wholeCenters(trajectory, k, ctx);
      // a frame of another topology (molecule count) or without a box: no continuity to keep
      if (!next.box || next.centers.length !== centers.length) {
        shifts = new Int32Array(next.centers.length);
      } else {
        // the jump between the frames in box vectors of the later one: rounded, it is a lattice move to undo
        // the later frame's box: this one going forward, the one before it (k + 1) going back
        const later = step > 0 ? next.box : (box ?? next.box);
        const inv = later.inverse;
        for (let m = 0, n = centers.length; m < n; m += 3) {
          // forward: n_k = n_(k-1) - round(inv (c_k - c_(k-1))); backward the same move undone
          const dx = (next.centers[m] - centers[m]) * step,
            dy = (next.centers[m + 1] - centers[m + 1]) * step,
            dz = (next.centers[m + 2] - centers[m + 2]) * step;
          const ra = Math.round(dx * inv[0] + dy * inv[3] + dz * inv[6]);
          const rb = Math.round(dx * inv[1] + dy * inv[4] + dz * inv[7]);
          const rc = Math.round(dx * inv[2] + dy * inv[5] + dz * inv[8]);
          shifts[m] -= step * ra;
          shifts[m + 1] -= step * rb;
          shifts[m + 2] -= step * rc;
        }
      }
      centers = next.centers;
      box = next.box;
      if (k % NoJumpEvery === 0 && !checkpoints.has(k)) checkpoints.set(k, Int32Array.from(shifts));
    }
    state.last = { index, shifts, centers, box };
    return shifts;
  }

  /** Each atom's molecule, by the forest's trees: parents are visited first, roots start molecules. */
  function withMolecules(plan: { parent: Int32Array; order?: Int32Array }): UnwrapPlan {
    const { parent, order } = plan;
    const molecule = new Int32Array(parent.length);
    let moleculeCount = 0;
    for (let k = 0, n = parent.length; k < n; k++) {
      const i = order ? order[k] : k;
      molecule[i] = parent[i] < 0 ? moleculeCount++ : molecule[parent[i]];
    }
    return { ...plan, molecule, moleculeCount };
  }

  /** The coordinates with every atom at the image nearest to its parent, in new arrays. */
  export function applyWhole(plan: UnwrapPlan, box: PeriodicBox, conformation: AtomicConformation): AtomicConformation {
    const x = Float32Array.from(conformation.x);
    const y = Float32Array.from(conformation.y);
    const z = Float32Array.from(conformation.z);
    const { parent, order } = plan;
    // the minimum image inlined (PeriodicBox.minimumImage): this runs once per atom per frame
    const [b0, b1, b2, b3, b4, b5, b6, b7, b8] = box.basis;
    const [i0, i1, i2, i3, i4, i5, i6, i7, i8] = box.inverse;
    for (let k = 0, n = parent.length; k < n; k++) {
      const i = order ? order[k] : k;
      const p = parent[i];
      if (p < 0) continue;
      const dx = x[i] - x[p],
        dy = y[i] - y[p],
        dz = z[i] - z[p];
      let fx = dx * i0 + dy * i3 + dz * i6;
      let fy = dx * i1 + dy * i4 + dz * i7;
      let fz = dx * i2 + dy * i5 + dz * i8;
      fx -= Math.round(fx);
      fy -= Math.round(fy);
      fz -= Math.round(fz);
      x[i] = x[p] + fx * b0 + fy * b3 + fz * b6;
      y[i] = y[p] + fx * b1 + fy * b4 + fz * b7;
      z[i] = z[p] + fx * b2 + fy * b5 + fz * b8;
    }
    return { ...conformation, id: UUID.create22(), x, y, z };
  }

  /** A BFS forest over the bonds: exact, whatever the frame's coordinates. */
  export function planFromBonds(bonds: IndexPairBonds, atomCount: number): UnwrapPlan {
    const { offset, b } = bonds.bonds;
    const { operatorA, operatorB } = bonds.bonds.edgeProps;
    const parent = new Int32Array(atomCount).fill(-1);
    const order = new Int32Array(atomCount);
    const seen = new Uint8Array(atomCount);
    let head = 0,
      tail = 0;
    for (let root = 0; root < atomCount; root++) {
      if (seen[root]) continue;
      seen[root] = 1;
      order[tail++] = root;
      while (head < tail) {
        const i = order[head++];
        for (let e = offset[i], _e = offset[i + 1]; e < _e; e++) {
          // a bond to another symmetry copy (mmCIF struct_conn) is not within the molecule
          if (bonds.hasOperators && operatorA[e] !== operatorB[e]) continue;
          const j = b[e];
          if (seen[j]) continue;
          seen[j] = 1;
          parent[j] = i;
          order[tail++] = j;
        }
      }
    }
    return withMolecules({ parent, order });
  }

  /**
   * From the hierarchy, not distances (the topology's own frame may be split): in a residue each atom's parent is the
   * atom before it; a polymer residue's first atom's parent is the first atom of the residue before it in the chain;
   * any other residue starts a molecule. Holds while these hops are under half a box (the minimum-image convention).
   */
  export function planFromHierarchy(hierarchy: AtomicHierarchy): UnwrapPlan {
    const { residueAtomSegments, chainAtomSegments, derived } = hierarchy;
    const { moleculeType } = derived.residue;
    const residueOffsets = residueAtomSegments.offsets;
    const parent = new Int32Array(hierarchy.atoms._rowCount).fill(-1);
    const residueOfChainStart = residueAtomSegments.index;
    for (let c = 0, _c = chainAtomSegments.count; c < _c; c++) {
      const rStart = residueOfChainStart[chainAtomSegments.offsets[c]];
      const rEnd = residueOfChainStart[chainAtomSegments.offsets[c + 1] - 1] + 1;
      for (let r = rStart; r < rEnd; r++) {
        const first = residueOffsets[r];
        if (r > rStart && isPolymer(moleculeType[r]) && isPolymer(moleculeType[r - 1])) {
          parent[first] = residueOffsets[r - 1];
        }
        for (let i = first + 1, _i = residueOffsets[r + 1]; i < _i; i++) parent[i] = i - 1;
      }
    }
    return withMolecules({ parent });
  }
}
