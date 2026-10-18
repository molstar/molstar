/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '@molstar/core/math/linear-algebra';
import type { ValueCell } from '@molstar/core/util';
import type { UniformValues } from './uniform.js';
import type { RenderableSchema, UniformSpec } from './renderable/schema.js';
import type { UniformKind } from './uniform.js';

/**
 * Depth is reversed on all backends: the near plane at depth 1, the far plane at depth 0 (the background, and the
 * cleared depth), so that floating-point depth keeps its precision far from the camera. The camera builds GL
 * projections (clip depth in `[-1, 1]`, near at -1); uniforms holding a projection are tagged (`UniformSpec.depthRange`)
 * and each backend maps them to its clip depth when it binds them, in double precision before they become `f32`:
 * WebGPU, and WebGL with `EXT_clip_control`, to `[0, 1]` reversed (depth equals clip z / w); WebGL without it to
 * `[-1, 1]` reversed (z negated; window depth is `(z / w + 1) / 2`, without the precision gain). Window depth is
 * `1 - d` of a GL projection's depth `d` on all.
 */
export type DepthRangeTag = 'projection' | 'inverse';

/** How a backend maps the clip coordinates of the camera's GL projections, see `DepthRangeTag` */
export interface DepthRangeConversion {
  /** a matrix that ends in the projection (`P`, `P * MV`) */
  readonly projection: (out: Mat4, m: Mat4) => Mat4;
  /** the inverse of such a matrix */
  readonly inverse: (out: Mat4, m: Mat4) => Mat4;
}

/** The conversion that multiplies clip coordinates by `clip` (column-major) */
export function createDepthRangeConversion(clip: Mat4): DepthRangeConversion {
  const inv = Mat4.invert(Mat4(), clip);
  return {
    projection: (out, m) => Mat4.mul(out, clip, m),
    inverse: (out, m) => Mat4.mul(out, m, inv),
  };
}

/** `[0, 1]` reversed: `z' = (w - z) / 2` */
export const ReversedZeroToOne = createDepthRangeConversion(
  Mat4.fromArray(Mat4(), [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -0.5, 0, 0, 0, 0.5, 1], 0),
);

/** `[-1, 1]` reversed: `z' = -z` */
export const ReversedMinusOneToOne = createDepthRangeConversion(
  Mat4.fromArray(Mat4(), [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1], 0),
);

/**
 * A cell with the converted value of `cell`, recomputed when its version changes; observers are those of `cell`
 * (`ValueCell.observe`)
 */
export function convertedCell(cell: ValueCell<Mat4>, convert: (out: Mat4, m: Mat4) => Mat4): ValueCell<Mat4> {
  const ref = { version: -1, value: Mat4() };
  return {
    get observers() {
      return cell.observers;
    },
    set observers(value) {
      cell.observers = value;
    },
    get ref() {
      if (cell.ref.version !== ref.version) {
        convert(ref.value, cell.ref.value);
        ref.version = cell.ref.version;
      }
      return ref;
    },
  } as unknown as ValueCell<Mat4>;
}

/** The cell to bind for uniform `name` of `schema`: tagged projections converted with `conversion` */
export function getDepthRangeCell<T>(
  schema: RenderableSchema,
  name: string,
  cell: ValueCell<T>,
  conversion: DepthRangeConversion | undefined,
): ValueCell<T> {
  const tag = conversion && (schema[name] as UniformSpec<UniformKind> | undefined)?.depthRange;
  if (!tag) return cell;
  return convertedCell(cell as unknown as ValueCell<Mat4>, conversion[tag]) as unknown as ValueCell<T>;
}

/** `values` with the tagged projections converted, see `getDepthRangeCell` */
export function convertDepthRange(
  schema: RenderableSchema,
  values: UniformValues,
  conversion: DepthRangeConversion | undefined,
): UniformValues {
  if (!conversion) return values;
  const converted: UniformValues = {};
  for (const k of Object.keys(values)) converted[k] = getDepthRangeCell(schema, k, values[k], conversion);
  return converted;
}
