/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';

// The members of a slot list's store (`ParticlePageStore`): a target its member while it has
// instances (its geometry, levels and bounds once), a freed member's index a new target's. What the store's owner keeps
// per member (a pool of values, `target/slots.ts`) is by these indices.

/** targets' ids under this: a member by a table, not a map */
const TableTargets = 1 << 22;

/** A target's member and its instances, by member index */
export class TargetMemberPool {
  /** a member's target, -1 a free member */
  targets = new Int32Array(64).fill(-1);
  /** a member's instances */
  counts = new Int32Array(64);
  /** per member its target's center (three a member) as a placement first gave it (`acquire`) */
  centers = new Float64Array(64 * 3);
  /** the members there are (used and free) */
  count = 0;
  /**
   * the members its owner holds room for (a pool): a target new to the store takes a member never used while there is
   * room, a free one's (of another target, the longest free) where there is none, one more past it where none is free
   */
  capacity = Infinity;

  /**
   * a target's member plus one (0 none) by target where its id is small, else a map: looked up a run (the yeast drawn
   * whole by its tiles: 2.29 M)
   */
  private ofTargetTable = new Int32Array(64);
  private readonly ofTargetMap = new Map<number, number>();

  private getOf(target: number): number | undefined {
    if (target >= 0 && target < TableTargets) {
      const m = target < this.ofTargetTable.length ? this.ofTargetTable[target] : 0;
      return m > 0 ? m - 1 : undefined;
    }
    return this.ofTargetMap.get(target);
  }

  private setOf(target: number, m: number) {
    if (target >= 0 && target < TableTargets) {
      if (target >= this.ofTargetTable.length) {
        const table = new Int32Array(Math.min(TableTargets, Math.max(target + 1, this.ofTargetTable.length * 2)));
        table.set(this.ofTargetTable);
        this.ofTargetTable = table;
      }
      this.ofTargetTable[target] = m + 1;
    } else this.ofTargetMap.set(target, m);
  }

  private deleteOf(target: number) {
    if (target >= 0 && target < TableTargets) {
      if (target < this.ofTargetTable.length) this.ofTargetTable[target] = 0;
    } else this.ofTargetMap.delete(target);
  }
  /** the target each member's center is of (-1 none): asked once a target's member, not once a particle's run */
  private centerTargets = new Int32Array(64).fill(-1);
  /** the free members, the longest free first (a `Set` keeps the order they were added in) */
  private readonly free = new Set<number>();
  /** a free member by the target it had: given to that target again before any other (its values that target's) */
  private readonly dormant = new Map<number, number>();
  /** a free member's last target */
  private readonly lastTargets = new Map<number, number>();

  /**
   * @param onGrow the members' arrays grew to `size`: the owner's per-member arrays grown alike
   * @param onAssign member `m` given to a target (its count 0): the owner's per-member values reset
   */
  constructor(
    private readonly onGrow?: (size: number) => void,
    private readonly onAssign?: (m: number) => void,
  ) {}

  /** `target`'s member, -1 where it has none */
  get(target: number) {
    return this.getOf(target) ?? -1;
  }

  /** `target`'s member, given one where it has none (its center `centerOf`'s, asked where the member had another) */
  acquire(target: number, centerOf: (target: number) => Vec3) {
    let m = this.getOf(target);
    if (m !== undefined) return m;
    // the member the target had where still free, else one never used while there is room, else the longest free
    m = this.dormant.get(target);
    if (m === undefined && this.count >= this.capacity && this.free.size > 0) {
      m = this.free.values().next().value!;
    }
    if (m !== undefined) {
      this.free.delete(m);
      this.dormant.delete(this.lastTargets.get(m)!);
      this.lastTargets.delete(m);
    } else {
      m = this.count++;
      if (m >= this.targets.length) this.grow(this.targets.length * 2);
    }
    this.targets[m] = target;
    this.counts[m] = 0;
    this.setOf(target, m);
    if (this.centerTargets[m] !== target) {
      const c = centerOf(target);
      this.centers[m * 3] = c[0];
      this.centers[m * 3 + 1] = c[1];
      this.centers[m * 3 + 2] = c[2];
      this.centerTargets[m] = target;
    }
    this.onAssign?.(m);
    return m;
  }

  /** Member `m` freed (its instances none): kept for its target's return, given to another where the pool is full */
  release(m: number) {
    const target = this.targets[m];
    this.deleteOf(target);
    this.targets[m] = -1;
    this.free.add(m);
    this.dormant.set(target, m);
    this.lastTargets.set(m, target);
  }

  private grow(size: number) {
    const targets = new Int32Array(size).fill(-1);
    targets.set(this.targets);
    this.targets = targets;
    const counts = new Int32Array(size);
    counts.set(this.counts);
    this.counts = counts;
    const centers = new Float64Array(size * 3);
    centers.set(this.centers);
    this.centers = centers;
    const centerTargets = new Int32Array(size).fill(-1);
    centerTargets.set(this.centerTargets);
    this.centerTargets = centerTargets;
    this.onGrow?.(size);
  }
}
