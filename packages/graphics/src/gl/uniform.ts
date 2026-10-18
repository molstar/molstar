/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Mat3, Mat4, Quat, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import type { ValueOf } from '@molstar/core/util/type-helpers';
import { deepClone } from '@molstar/core/util/object';
import type { UniformBlockWriter } from './uniform-layout.js';

export type UniformKindValue = {
  b: boolean;
  'b[]': boolean[];
  f: number;
  'f[]': number[];
  i: number;
  'i[]': number[];
  v2: Vec2;
  'v2[]': number[];
  v3: Vec3;
  'v3[]': number[];
  v4: Vec4;
  'v4[]': number[];
  q: Quat;
  'q[]': number[];
  iv2: Vec2;
  'iv2[]': number[];
  iv3: Vec3;
  'iv3[]': number[];
  iv4: Vec4;
  'iv4[]': number[];
  m3: Mat3;
  'm3[]': number[];
  m4: Mat4;
  'm4[]': number[];
  t: number;
  't[]': number[];
};
export type UniformKind = keyof UniformKindValue;
export type UniformType = ValueOf<UniformKindValue>;

export type UniformValues = { [k: string]: ValueCell<UniformType> };
export type UniformsList = [string, ValueCell<UniformType>][];

export function isArrayUniform(kind: UniformKind) {
  return kind.endsWith('[]');
}

export function isUniformValueScalar(kind: UniformKind): boolean {
  switch (kind) {
    case 'f':
    case 'i':
    case 'b':
      return true;
    default:
      return false;
  }
}

export function cloneUniformValues(uniformValues: UniformValues): UniformValues {
  const clonedValues: UniformValues = {};
  Object.keys(uniformValues).forEach((k) => {
    clonedValues[k] = ValueCell.create(deepClone(uniformValues[k].ref.value));
  });
  return clonedValues;
}

/** Instance of a uniform block layout, see `getUniformBlockLayout` */
export interface UniformBuffer {
  readonly id: number;

  getByteCount: () => number;
  /**
   * Create a writer that copies the given uniform values into the CPU-side
   * data. Each call to the writer only copies values whose version changed.
   * Cell instances are resolved once, keep using the same writer for them.
   */
  createWriter: (uniformValues: UniformsList) => UniformBlockWriter;
  /** upload CPU-side changes */
  sync: () => void;
  /** bind the buffer at its block's fixed binding index */
  bind: () => void;

  destroy: () => void;
}
