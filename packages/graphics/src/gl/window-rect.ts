/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec2, Vec4 } from '@molstar/core/math/linear-algebra';
import type { ValueCell } from '@molstar/core/util';
import type { RenderableSchema, UniformSpec } from './renderable/schema.js';
import type { UniformKind, UniformValues } from './uniform.js';

/**
 * Window rectangles have a top-left origin (y down) in the neutral API, as in the browser; uniforms holding one are
 * tagged (`UniformSpec.windowRect`) so that a backend with a bottom-left origin (WebGL, whose fragment coordinates and
 * rendered textures are bottom-up) flips them where it binds them:
 * - `bounds`: `vec4(x0, y0, x1, y1)` in uv of the target
 * - `rect`: `vec4(x, y, width, height)` in uv of the target
 * - `pixels`: `vec4(x, y, width, height)` in pixels of the target of the render pass
 * - `offset`: `vec2(x, y)`, a displacement in uv
 */
export type WindowRectTag = 'bounds' | 'rect' | 'pixels' | 'offset';

/** `v` flipped to a bottom-left origin in a target of `height` pixels (`pixels`) or uv (the others) */
export function flipWindowRect(out: Vec4, v: Vec4, tag: WindowRectTag, height: number): Vec4 {
  switch (tag) {
    case 'offset':
      out[0] = v[0];
      out[1] = -v[1];
      return out;
    case 'bounds':
      return Vec4.set(out, v[0], 1 - v[3], v[2], 1 - v[1]);
    case 'rect':
      return Vec4.set(out, v[0], 1 - v[1] - v[3], v[2], v[3]);
    case 'pixels':
      return Vec4.set(out, v[0], height - v[1] - v[3], v[2], v[3]);
  }
}

/**
 * A cell with the flipped value of `cell`, recomputed when its version changes or, for `pixels`, the height of the
 * render pass target (`getTargetHeight`); observers are those of `cell` (`ValueCell.observe`), and `pixels` cells are
 * volatile (the target height changes without a notification)
 */
function flippedCell(cell: ValueCell<Vec4>, tag: WindowRectTag, getTargetHeight: () => number): ValueCell<Vec4> {
  const ref = { version: 0, value: tag === 'offset' ? (Vec2() as unknown as Vec4) : Vec4() };
  let version = -1;
  let height = -1;
  return {
    get observers() {
      return cell.observers;
    },
    set observers(value) {
      cell.observers = value;
    },
    volatile: tag === 'pixels',
    get ref() {
      const h = tag === 'pixels' ? getTargetHeight() : 0;
      if (cell.ref.version !== version || h !== height) {
        flipWindowRect(ref.value, cell.ref.value, tag, h);
        version = cell.ref.version;
        height = h;
        ref.version += 1;
      }
      return ref;
    },
  } as unknown as ValueCell<Vec4>;
}

/** The cell to bind for uniform `name` of `schema`: tagged window rectangles flipped */
export function getWindowRectCell<T>(
  schema: RenderableSchema,
  name: string,
  cell: ValueCell<T>,
  getTargetHeight: (() => number) | undefined,
): ValueCell<T> {
  const tag = getTargetHeight && (schema[name] as UniformSpec<UniformKind> | undefined)?.windowRect;
  if (!tag) return cell;
  return flippedCell(cell as unknown as ValueCell<Vec4>, tag, getTargetHeight) as unknown as ValueCell<T>;
}

/** `values` with the tagged window rectangles flipped, see `getWindowRectCell` */
export function flipWindowRects(
  schema: RenderableSchema,
  values: UniformValues,
  getTargetHeight: (() => number) | undefined,
): UniformValues {
  if (!getTargetHeight) return values;
  const flipped: UniformValues = {};
  for (const k of Object.keys(values)) flipped[k] = getWindowRectCell(schema, k, values[k], getTargetHeight);
  return flipped;
}
