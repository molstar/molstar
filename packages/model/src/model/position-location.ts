/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';

/** A position location used by geometry that samples a position and normal. */
export interface PositionLocation {
  readonly kind: 'position-location';
  readonly position: Vec3;
  /** Normal vector at the position (used for surface coloring) */
  readonly normal: Vec3;
}

export function PositionLocation(position?: Vec3, normal?: Vec3): PositionLocation {
  return {
    kind: 'position-location',
    position: position ? Vec3.clone(position) : Vec3(),
    normal: normal ? Vec3.clone(normal) : Vec3(),
  };
}

export function isPositionLocation(x: any): x is PositionLocation {
  return !!x && x.kind === 'position-location';
}
