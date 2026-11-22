/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

/**
 * antimatter15's gaussian splats (`.splat`): no header, a row of `SplatRowBytes` a splat, little-endian. Its bytes as
 * given (no copy), read through `view` at the `SplatRow` offsets
 */
export interface SplatFile {
  readonly count: number;
  readonly data: ReadonlyUint8Array;
  readonly view: DataView;
}

export const SplatRowBytes = 32;

/** Byte offsets within a row */
export const SplatRow = {
  /** x, y, z: `float32` */
  position: 0,
  /** the scales along the splat's axes (linear): `float32` */
  scale: 12,
  /** r, g, b (the base color) and the opacity (sigmoid's): bytes */
  color: 24,
  /** the rotation `w, x, y, z` as bytes `q * 128 + 128` */
  rotation: 28,
} as const;
