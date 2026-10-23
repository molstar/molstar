/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PeriodicBox } from '@molstar/core/math/geometry/periodic-box';
import { UUID } from '@molstar/core/util/uuid';
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';
import { Model } from './model.js';
import type { AtomicConformation, AtomicHierarchy, AtomicRanges } from './properties/atomic.js';
import { getAtomicRanges } from './properties/utils/atomic-ranges.js';
import { isPolymer } from './types.js';

// Molecules split by a simulation's periodic box made whole again (.v6/periodic-trajectories/plan.md §3): every atom
// moved by a lattice vector to the image nearest to its parent in a spanning forest of the topology, built once.

export { ModelPeriodicity };

/**
 * How to make a topology's molecules whole: each atom's parent, the atom it is imaged next to (-1: a root, it stays
 * where it is), and the order to visit them in, parents first.
 */
interface UnwrapPlan {
  readonly parent: Int32Array;
  /** Atoms in visiting order; absent: in index order (every parent comes before its children) */
  readonly order?: Int32Array;
  /** The topology's ranges computed on a whole frame (polymer gaps): shared by its frames once made whole */
  ranges?: AtomicRanges;
}

namespace ModelPeriodicity {
  const plans = new WeakMap<object, UnwrapPlan>();
  const wholeFrames = new WeakMap<Model, Model>();

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
   * The frame with its molecules whole: a new model (the same one for the same frame), or the frame itself when it
   * has no simulation box (`isApplicable`).
   */
  export function makeWhole(model: Model): Model {
    const box = boxOf(model);
    if (!box) return model;
    const cached = wholeFrames.get(model);
    if (cached) return cached;

    const plan = getPlan(model);
    const conformation = applyWhole(plan, box, model.atomicConformation);
    // polymer gaps were found on the topology's coordinates, which may have been split: found again on whole ones
    if (!plan.ranges) {
      plan.ranges = getAtomicRanges(model.atomicHierarchy, model.entities, conformation, model.sequence);
    }
    const whole = Model.withAtomicConformation(model, conformation, plan.ranges);
    wholeFrames.set(model, whole);
    return whole;
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
    return { parent, order };
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
    return { parent };
  }
}
