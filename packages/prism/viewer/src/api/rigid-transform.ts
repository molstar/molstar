/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { TransformParam } from '@molstar/plugin/state/transforms/helpers';
import { type ApiFailure, ApiResult } from './area.js';
import { OptionalParam } from './params.js';

// A move or turn of a whole object by hand (`structure.transform`, `particles.transform`): the same arguments, a
// matrix or parts (a translation, an angle about an axis), from where the object was loaded or on top of where it is.

/** The params of a move: `what` names the object in the center's options ("structure") */
export function rigidTransformParams(what: string) {
  return {
    translation: OptionalParam(PD.Vec3(Vec3.create(0, 0, 0), {}, { description: 'Translation in Å: [x, y, z]' })),
    axis: OptionalParam(
      PD.Vec3(
        Vec3.create(0, 0, 1),
        {},
        { description: 'Rotation axis: [x, y, z], with an angle ([0, 0, 1] when absent)' },
      ),
    ),
    angle: OptionalParam(PD.Numeric(0, { min: -360, max: 360, step: 1 }, { description: 'Rotation angle in degrees' })),
    center: PD.Select<'centroid' | 'origin'>(
      'centroid',
      [
        ['centroid', `The ${what}’s center`],
        ['origin', 'The origin'],
      ],
      { description: `What the rotation turns about: the ${what}’s center or the origin` },
    ),
    matrix: OptionalParam(
      PD.Mat4(Mat4.identity(), {
        description: 'A 4×4 matrix (column-major, 16 numbers) instead of the axis, angle and translation',
      }),
    ),
    compose: PD.Boolean(false, { description: 'Apply on top of the current transform instead of replacing it' }),
  };
}

/** What a move's arguments ask for (api.md §1: a matrix or parts, not both, not none); `reset` names the reset */
export function wantedTransform(v: { center: 'centroid' | 'origin' }, reset: string): TransformParam | ApiFailure {
  const {
    translation: t,
    axis: a,
    angle: g,
    matrix,
  } = v as {
    translation?: Vec3;
    axis?: Vec3;
    angle?: number;
    matrix?: Mat4;
  };
  const byMatrix = matrix !== undefined;
  const byParts = t !== undefined || g !== undefined || a !== undefined;
  if (byMatrix && byParts) {
    return ApiResult.fail(
      'invalid',
      'Give a matrix, or an axis, angle and translation, not both.',
      'Leave the matrix out to move by parts.',
    );
  }
  if (!byMatrix && !byParts) {
    return ApiResult.fail(
      'invalid',
      'Move how?',
      `Give a translation, an angle (about an axis), or a matrix; ${reset} puts it back.`,
    );
  }
  if (a !== undefined && g === undefined) {
    return ApiResult.fail('invalid', 'An axis needs an angle.', 'Give the angle in degrees with the axis.');
  }
  const translation = t ?? Vec3.zero();
  const axis = a ?? Vec3.create(0, 0, 1);
  const angle = g ?? 0;
  if (!byMatrix && angle !== 0 && Vec3.magnitude(axis) < 1e-6) {
    return ApiResult.fail('invalid', 'The rotation axis is zero.', 'Give an axis such as [0, 0, 1].');
  }
  return byMatrix
    ? { name: 'matrix', params: { data: Mat4.clone(matrix!), transpose: false } }
    : {
        name: 'components',
        params: {
          translation: Vec3.clone(translation),
          axis: Vec3.normalize(Vec3(), Vec3.magnitude(axis) < 1e-6 ? Vec3.create(0, 0, 1) : axis),
          angle,
          rotationCenter:
            v.center === 'centroid'
              ? { name: 'centroid', params: {} }
              : { name: 'point', params: { point: Vec3.zero() } },
        },
      };
}

/** The transform that puts an object back where it was loaded (a decorator kept: removing it takes its subtree) */
export const IdentityTransform: TransformParam = {
  name: 'components',
  params: {
    translation: Vec3.zero(),
    axis: Vec3.create(0, 0, 1),
    angle: 0,
    rotationCenter: { name: 'centroid', params: {} },
  },
};
