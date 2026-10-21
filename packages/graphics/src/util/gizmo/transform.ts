/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Color } from '@molstar/core/util/color';
import type { Mesh } from '../../geo/geometry/mesh/mesh.js';
import { MeshBuilder } from '../../geo/geometry/mesh/mesh-builder.js';
import { Box } from '../../geo/primitive/box.js';
import { AxisColors, type Handle, HandleColors, type HandleGroup } from './handle.js';
import { arrow, knob, ring } from './primitives.js';

// The standard transform handles, as in 3D tools: an arrow along each axis (move along it), a square in each pair of
// axes' plane (move in it), a ring around each axis (turn about it), the center (move in the screen's plane); a knob
// past each arrow's tip scales along its axis. In the item's frame, `size` the arrows' length; X red, Y green, Z blue.

export interface TransformHandlesOptions {
  /** Which parts (all but `scale` by default) */
  readonly move?: boolean;
  readonly planes?: boolean;
  readonly rotate?: boolean;
  readonly scale?: boolean;
}

const Axes = [Vec3.create(1, 0, 0), Vec3.create(0, 1, 0), Vec3.create(0, 0, 1)] as const;
const Names = ['X', 'Y', 'Z'] as const;

/** The transform handles' mesh (one group a handle) and what each does */
export function transformHandles(
  size: number,
  options: TransformHandlesOptions = {},
): { handles: Mesh; groups: HandleGroup[] } {
  const { move = true, planes = true, rotate = true, scale = false } = options;
  const state = MeshBuilder.createState(4096, 2048);
  const groups: HandleGroup[] = [];
  const group = (handle: Handle, label: string, color: Color) => {
    state.currentGroup = groups.length;
    groups.push({ handle, label, color });
  };
  const shaft = size * 0.02;
  if (move) {
    group({ kind: 'move' }, 'Drag to move it', HandleColors.move);
    knob(state, Vec3(), size * 0.07);
    for (let i = 0; i < 3; i++) {
      group({ kind: 'slide', axis: Axes[i] }, `Drag to move it along ${Names[i]}`, AxisColors[i]);
      arrow(state, Vec3.scale(Vec3(), Axes[i], size * 0.1), Vec3.scale(Vec3(), Axes[i], size), shaft);
    }
  }
  if (planes) {
    // the plane across axis i, between the other two: a thin square near the origin
    for (let i = 0; i < 3; i++) {
      const [b, c] = [0, 1, 2].filter((x) => x !== i);
      const at = Vec3();
      at[b] = at[c] = size * 0.3;
      const extent = Vec3.create(size * 0.2, size * 0.2, size * 0.2);
      extent[i] = size * 0.01;
      group({ kind: 'plane', normal: Axes[i] }, `Drag to move it in ${Names[b]}${Names[c]}`, AxisColors[i]);
      MeshBuilder.addPrimitive(
        state,
        Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), at), Mat4.fromScaling(Mat4(), extent)),
        Box(),
      );
    }
  }
  if (rotate) {
    // a ring around axis i: around +Y as built, turned onto X (about Z) and Z (about X)
    const onto = [Mat4.rotZ90 as Mat4, Mat4.identity(), Mat4.rotX90 as Mat4];
    for (let i = 0; i < 3; i++) {
      group({ kind: 'rotate', axis: Axes[i] }, `Drag to turn it about ${Names[i]}`, AxisColors[i]);
      ring(state, size * 0.8, shaft * 0.8, 0, Mat4.mul(Mat4(), onto[i], Mat4.rotX90 as Mat4));
    }
  }
  if (scale) {
    for (let i = 0; i < 3; i++) {
      group({ kind: 'scale', axes: [i as 0 | 1 | 2] }, `Drag to size it along ${Names[i]}`, AxisColors[i]);
      knob(state, Vec3.scale(Vec3(), Axes[i], size * 1.12), size * 0.04);
    }
  }
  return { handles: MeshBuilder.getMesh(state), groups };
}
