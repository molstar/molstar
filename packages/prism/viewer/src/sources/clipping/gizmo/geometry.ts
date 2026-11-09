/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Color } from '@molstar/core/util/color';
import { clipObjectFrame, clipWallsMesh } from '@molstar/graphics/geo/clip-volumes';
import type { Lines } from '@molstar/graphics/geo/geometry/lines/lines';
import { LinesBuilder } from '@molstar/graphics/geo/geometry/lines/lines-builder';
import type { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import type { Points } from '@molstar/graphics/geo/geometry/points/points';
import { PointsBuilder } from '@molstar/graphics/geo/geometry/points/points-builder';
import { type Handle, HandleColors, type HandleGroup } from '@molstar/graphics/util/gizmo/handle';
import type { GizmoPart } from '@molstar/graphics/util/gizmo/controller';
import { ellipseLines, ringLines } from '@molstar/graphics/util/gizmo/primitives';
import { type TransformHandlesOptions, transformHandles } from '@molstar/graphics/util/gizmo/transform';
import { type CutObject, coneAngle, LocalUp } from '../cuts.js';

// A cut's gizmo (clipping.md §4): its body (the clip volume's walls, `graphics/geo/clip-volumes`) and its handles, one
// group each saying what a drag does (the gizmo library's `Handle`). Built in the cut's frame (a plane's from where it
// crosses what it cuts), so a move or a turn only moves the frame. In two parts:
// - the shape, at the scene's scale (sized by what it cuts: `reach`, the radius of their bounds): the clip object and
//   the handles that size it, a box's edges, corners and face knobs, an ellipsoid's (the sphere's) sections, tips and
//   octant knobs, a cylinder's (elliptic) rims, side knobs and cap, a cone's base rim; edges and rims as lines (as
//   wide on screen at any zoom);
// - the move and turn handles, at a fixed size on screen (`MoveScreenSize`, as the move gizmo): the gizmo
//   library's transform handles (Blender's move and rotate gizmo) about the frame's origin, by cut (`moveOptions`).
// The binding gives them to the controller as parts (`GizmoLayout.parts`); the move handles' groups first.

/** The cone's angle handle (a `custom` handle, the binding's own gesture): its base rim widened or narrowed */
export const ConeAngleHandle = 'cone-angle';

/** The move handles' size on screen: CSS pixels a unit (an arrow is 1 long, the rings' radius 0.65) */
export const MoveScreenSize = 90;

export interface GizmoGeometry {
  readonly body: Mesh;
  /** The shape's edges and rims as lines, a set width on screen */
  readonly handleLines?: Lines;
  /** The shape's knobs as points, a set size on screen */
  readonly handlePoints?: Points;
  readonly groups: readonly HandleGroup[];
  /** The move and turn handles' parts (`transformHandles`) at `MoveScreenSize` */
  readonly moveParts: readonly GizmoPart[];
  /** The cut's frame in the world */
  readonly transform: Mat4;
}

/** The meshes of a cut's gizmo over what it cuts */
export function gizmoGeometry(object: CutObject, bounds: Sphere3D): GizmoGeometry {
  const reach = Math.max(bounds.radius, 1);
  const move = transformHandles(1, moveOptions(object));
  const groups: HandleGroup[] = [...move.groups];
  const shape = buildShape(object, reach, groups);
  return {
    // the walls as the clip volumes shape draws them
    body: clipWallsMesh(object, reach),
    ...(shape.lines ? { handleLines: shape.lines } : {}),
    ...(shape.points ? { handlePoints: shape.points } : {}),
    groups,
    moveParts: move.parts.map((p) => ({ ...p, screenSize: MoveScreenSize })),
    transform: clipObjectFrame(object, bounds),
  };
}

/**
 * A cut's move and turn handles (`transformHandles`): what changes it. A plane's center slides along its normal and
 * its one arrow points to the side it cuts away (moving it within itself changes nothing: no squares), its rings turn
 * its normal; a round cylinder or a cone no ring about its own axis (round about it: it changes nothing), an elliptic
 * cylinder one (its radii, the scale's x and z, turn with it)
 */
/** Whether a cylinder of this scale is round: its radii (x, z) equal within a thousandth */
function isRound(scale: Vec3) {
  return Math.abs(scale[0] - scale[2]) <= 1e-3 * Math.max(scale[0], scale[2]);
}

function moveOptions(object: CutObject): TransformHandlesOptions {
  switch (object.type) {
    case 'plane':
      return {
        centerAxis: LocalUp,
        arrows: [Vec3.scale(Vec3(), LocalUp, object.invert ? -1 : 1)],
        planes: false,
        rings: [0, 2],
      };
    case 'cylinder':
      return isRound(object.scale) ? { rings: [0, 2] } : {};
    // its axis its opening, local −Z
    case 'infiniteCone':
      return { rings: [0, 1] };
    default:
      return {};
  }
}

/**
 * A cut's shape handles: its edges and rims as lines and its knobs as points (at their places in the world, a set size
 * on screen), their groups into `groups` after the move handles'
 */
function buildShape(object: CutObject, reach: number, groups: HandleGroup[]) {
  const lines = LinesBuilder.create(512, 256);
  const points = PointsBuilder.create(16, 16);
  let current = 0;
  const group = (handle: Handle, label: string, color: Color) => {
    current = groups.length;
    groups.push({ handle, label, color });
  };
  // the shape's parts, drawn only
  const edgeLine = (from: Vec3, to: Vec3) => lines.add(from[0], from[1], from[2], to[0], to[1], to[2], current);
  // a rim of radius `r` around the frame's up (+Y), placed by `t`
  const rim = (r: number, t?: Mat4) => ringLines(lines, r, current, t);
  const knobAt = (at: Vec3) => points.add(at[0], at[1], at[2], current);
  switch (object.type) {
    case 'plane':
      // its move and turn handles alone
      break;
    case 'sphere': {
      // an ellipsoid (the shader's sphere, its half axes the scale's halves): the box's handles on its surface
      const half = [object.scale[0] / 2, object.scale[1] / 2, object.scale[2] / 2];
      // its tips: a knob at each moves that tip along its axis (the opposite one stays), as a box's face knobs
      for (const a of [0, 1, 2] as const) {
        for (const sign of [-1, 1] as const) {
          const at = Vec3();
          at[a] = sign * half[a];
          group({ kind: 'faces', faces: [{ axis: a, sign }] }, 'Drag to move this side', HandleColors.outline);
          knobAt(at);
        }
      }
      // its sections through two axes, the edges joining the tips: each sizes those two, as a box's edges
      const sections: [0 | 1 | 2, 0 | 1 | 2, Mat4 | undefined][] = [
        [0, 2, undefined],
        [0, 1, Mat4.rotX90 as Mat4],
        [1, 2, Mat4.rotZ90 as Mat4],
      ];
      for (const [b, c, t] of sections) {
        group({ kind: 'scale', axes: [b, c] }, 'Drag to size it across this section', HandleColors.outline);
        ellipseLines(lines, half[b], half[c], current, t);
      }
      // a knob at each octant's middle on its surface (a box corner's place) sizes it whole
      group({ kind: 'scale' }, 'Drag to size it whole', HandleColors.outline);
      for (const x of [-1, 1]) {
        for (const y of [-1, 1]) {
          for (const z of [-1, 1]) {
            knobAt(Vec3.scale(Vec3(), Vec3.create(x * half[0], y * half[1], z * half[2]), 1 / Math.sqrt(3)));
          }
        }
      }
      break;
    }
    case 'cube': {
      const [sx, sy, sz] = object.scale;
      const half = [sx / 2, sy / 2, sz / 2];
      // its faces: a knob at the center of each moves that face along its normal (the opposite one stays)
      for (const a of [0, 1, 2] as const) {
        for (const sign of [-1, 1] as const) {
          const at = Vec3();
          at[a] = sign * half[a];
          group({ kind: 'faces', faces: [{ axis: a, sign }] }, 'Drag to move this side', HandleColors.outline);
          knobAt(at);
        }
      }
      // its edges, the four along each axis one handle (lit together): the section across that axis scaled
      for (const a of [0, 1, 2] as const) {
        const [b, c] = ([0, 1, 2] as const).filter((x) => x !== a);
        group({ kind: 'scale', axes: [b, c] }, 'Drag to size it across these edges', HandleColors.outline);
        for (const sb of [-1, 1]) {
          for (const sc of [-1, 1]) {
            const from = Vec3();
            const to = Vec3();
            from[a] = -half[a];
            to[a] = half[a];
            from[b] = to[b] = sb * half[b];
            from[c] = to[c] = sc * half[c];
            edgeLine(from, to);
          }
        }
      }
      // its corners scale it whole
      group({ kind: 'scale' }, 'Drag to size it whole', HandleColors.outline);
      for (const x of [-1, 1]) {
        for (const y of [-1, 1]) {
          for (const z of [-1, 1]) knobAt(Vec3.create(x * half[0], y * half[1], z * half[2]));
        }
      }
      break;
    }
    case 'cylinder': {
      // elliptic (the shader's: its radii the scale's x and z halves)
      const [a, h, c] = [object.scale[0] / 2, object.scale[1], object.scale[2] / 2];
      // its rims size both radii (about their own centers), in proportion
      for (const y of [-h / 2, h / 2]) {
        group(
          { kind: 'scale', axes: [0, 2], pivot: Vec3.create(0, y, 0) },
          'Drag to size its radii',
          HandleColors.outline,
        );
        ellipseLines(lines, a, c, current, Mat4.fromTranslation(Mat4(), Vec3.create(0, y, 0)));
      }
      // a knob on each side halfway up moves that side across the axis (the opposite one stays), as a box's face knobs
      for (const [axis, r] of [
        [0, a],
        [2, c],
      ] as const) {
        for (const sign of [-1, 1] as const) {
          const at = Vec3();
          at[axis] = sign * r;
          group({ kind: 'faces', faces: [{ axis, sign }] }, 'Drag to move this side', HandleColors.outline);
          knobAt(at);
        }
      }
      // a knob at each cap's center sizes its height (the other cap stays)
      for (const sign of [-1, 1] as const) {
        group({ kind: 'faces', faces: [{ axis: 1, sign }] }, 'Drag to move this end', HandleColors.outline);
        knobAt(Vec3.create(0, (sign * h) / 2, 0));
      }
      break;
    }
    case 'infiniteCone': {
      // drawn over what it cuts: from its apex as far as their reach, along its opening (local −Z)
      const length = reach;
      const radius = length * Math.tan((coneAngle(object.scale) * Math.PI) / 180);
      // its base rim widens or narrows it (about the base's center), around its opening
      group(
        { kind: 'custom', id: ConeAngleHandle, pivot: Vec3.create(0, 0, -length) },
        'Drag to widen or narrow it',
        HandleColors.outline,
      );
      rim(radius, Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, -length)), Mat4.rotX90 as Mat4));
      break;
    }
  }
  const builtLines = lines.getLines();
  const builtPoints = points.getPoints();
  return {
    lines: builtLines.lineCount > 0 ? builtLines : undefined,
    points: builtPoints.pointCount > 0 ? builtPoints : undefined,
  };
}
