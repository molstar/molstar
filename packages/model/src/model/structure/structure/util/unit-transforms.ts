/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Structure, Unit } from '../../../structure.js';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { IntMap } from '@molstar/core/data/int';
import { fillIdentityTransform } from '@molstar/core/math/linear-algebra/3d/transform-array';

const tmpMat = Mat4();
const identity = Mat4.identity();

export class StructureUnitTransforms {
  private unitTransforms: Float32Array;
  private groupUnitTransforms: Float32Array[] = [];
  /** maps unit.id to offset of transform in unitTransforms */
  private unitOffsetMap = IntMap.Mutable<number>();
  private groupIndexMap = new Map<Unit.SymmetryGroup, number>();
  private size: number;

  private _isIdentity: boolean | undefined = undefined;

  version = 0;

  constructor(readonly structure: Structure) {
    this.unitTransforms = new Float32Array(structure.units.length * 16);
    this.size = structure.units.length;
    this.reset(); // to set to identity
    let groupOffset = 0;
    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
      const g = structure.unitSymmetryGroups[i];
      this.groupIndexMap.set(g, i);
      const groupTransforms = this.unitTransforms.subarray(groupOffset, groupOffset + g.units.length * 16);
      this.groupUnitTransforms.push(groupTransforms);
      for (let j = 0, jl = g.units.length; j < jl; ++j) {
        this.unitOffsetMap.set(g.units[j].id, groupOffset + j * 16);
      }
      groupOffset += g.units.length * 16;
    }
  }

  reset() {
    this.version = 0;
    fillIdentityTransform(this.unitTransforms, this.size);
    this._isIdentity = true;
  }

  get isIdentity() {
    if (this._isIdentity === undefined) {
      this._isIdentity = true;
      for (let i = 0, il = this.size * 16; i < il; i += 16) {
        Mat4.fromArray(tmpMat, this.unitTransforms, i);
        if (!Mat4.isIdentity(tmpMat)) {
          this._isIdentity = false;
          break;
        }
      }
    }
    return this._isIdentity;
  }

  setTransform(matrix: Mat4, unit: Unit) {
    this.version = (this.version + 1) % 0x7fffffff;
    Mat4.toArray(matrix, this.unitTransforms, this.unitOffsetMap.get(unit.id));
    this._isIdentity = undefined;
  }

  getTransform(out: Mat4, unit: Unit) {
    return Mat4.fromArray(out, this.unitTransforms, this.unitOffsetMap.get(unit.id));
  }

  /** per group not of `structure.unitSymmetryGroups` its units' transforms, gathered by unit id at each call */
  private otherGroupTransforms = new WeakMap<Unit.SymmetryGroup, Float32Array>();

  /**
   * The transforms of `group`'s units, in its order: a view of the structure's own groups' (`unitSymmetryGroups`), and
   * for any other group (a representation's, e.g. `computeChildAwareTransformGroups`: new group objects) gathered by
   * unit id; a unit not of the structure identity
   */
  getSymmetryGroupTransforms(group: Unit.SymmetryGroup): Float32Array {
    const index = this.groupIndexMap.get(group);
    if (index !== undefined) return this.groupUnitTransforms[index];
    const { units } = group;
    let out = this.otherGroupTransforms.get(group);
    if (!out) this.otherGroupTransforms.set(group, (out = new Float32Array(units.length * 16)));
    for (let j = 0, jl = units.length; j < jl; ++j) {
      const offset = this.unitOffsetMap.get(units[j].id);
      if (offset === undefined) Mat4.toArray(identity, out, j * 16);
      else out.set(this.unitTransforms.subarray(offset, offset + 16), j * 16);
    }
    return out;
  }
}
