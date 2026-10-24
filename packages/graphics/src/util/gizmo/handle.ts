/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import type { Vec3 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';

// A gizmo's handles: what dragging each does (`Handle`), how it reads (its label, its color). A gizmo's handle mesh has
// one group a handle; the controller finds a hovered handle by its group (`HandleGroup`, the group's id its index).

/** A face of a box (or a cylinder's end): its local axis and which way it faces */
export interface Face {
  readonly axis: 0 | 1 | 2;
  readonly sign: 1 | -1;
}

/**
 * What dragging a handle does to its item's pose (`GizmoPose`):
 * - `move`: in the screen's plane;
 * - `slide`: along `axis` (local), as the screen shows it;
 * - `plane`: in the plane across `normal` (local);
 * - `turn`: about the screen's axes;
 * - `rotate`: about `axis` (local), the point grabbed on the ring following the cursor (the cursor taken onto the ring's
 *   plane; along the ring's tangent there where the plane is edge-on), so either side of a ring turns the cursor's way;
 * - `scale`: `axes` (all when absent) by the cursor's distance from `pivot` (local; the item's origin when absent) over
 *   the grab's: the part held stays under the cursor whichever side it is on;
 * - `faces`: faces dragged along their axes, the opposite ones staying (a box's face, a cylinder's end);
 * - `custom`: the binding's own (a cone's angle), its `pivot` for the grab as `scale`'s.
 */
export type Handle =
  | { readonly kind: 'move' }
  | { readonly kind: 'slide'; readonly axis: Vec3 }
  | { readonly kind: 'plane'; readonly normal: Vec3 }
  | { readonly kind: 'turn' }
  | { readonly kind: 'rotate'; readonly axis: Vec3 }
  | { readonly kind: 'scale'; readonly axes?: readonly (0 | 1 | 2)[]; readonly pivot?: Vec3 }
  | { readonly kind: 'faces'; readonly faces: readonly Face[] }
  | { readonly kind: 'custom'; readonly id: string; readonly pivot?: Vec3 };

/** A handle, how it reads (a tooltip), its color */
export interface HandleGroup {
  readonly handle: Handle;
  readonly label: string;
  readonly color: Color;
}

/** The colors by what a handle does: move, slide, turn; sizing parts (outlines, knobs on faces) one light accent */
export const HandleColors = {
  move: Color(0xf5a623),
  slide: Color(0x4a90e2),
  turn: Color(0x50b860),
  outline: Color(0x9db8e0),
} as const;

/** The axes' colors, as everywhere in 3D tools: X red, Y green, Z blue */
export const AxisColors = [Color(0xe05a5a), Color(0x5cbf5c), Color(0x4a90e2)] as const;
