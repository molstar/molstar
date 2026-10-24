/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import type { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Color } from '@molstar/core/util/color';
import { clipObjectFrame, clipWallsMesh } from '@molstar/graphics/geo/clip-volumes';
import type { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { type Handle, HandleColors, type HandleGroup } from '@molstar/graphics/util/gizmo/handle';
import { arrow, knob, ring, tube } from '@molstar/graphics/util/gizmo/primitives';
import { ConeOpening, type CutObject, coneAngle, LocalUp } from '../cuts.js';

// A cut's gizmo (clipping.md §4): its body (the clip volume's walls, `graphics/geo/clip-volumes`) and its handles, one
// group each saying what a drag does (the gizmo library's `Handle`). Built in the cut's frame (a plane's from where it
// crosses what it cuts), so a move or a turn only moves the frame. Sized by what it cuts (`reach`, the radius of their
// bounds): knobs, tubes and rings stay in proportion to the scene. A closed shape's outline is its handles: a box's
// edges and corners (and a knob a face), a cylinder's rims, a cone's base rim, a sphere's ring (facing the camera).

/** The cone's angle handle (a `custom` handle, the binding's own gesture): its base rim widened or narrowed */
export const ConeAngleHandle = 'cone-angle';

export interface GizmoGeometry {
  readonly body: Mesh;
  readonly handles: Mesh;
  readonly groups: readonly HandleGroup[];
  /** The cut's frame in the world */
  readonly transform: Mat4;
  /** The frame turned to face the camera (a sphere: its ring always round) */
  readonly facesCamera: boolean;
}

/** The meshes of a cut's gizmo over what it cuts */
export function gizmoGeometry(object: CutObject, bounds: Sphere3D): GizmoGeometry {
  const reach = Math.max(bounds.radius, 1);
  const knobR = reach * 0.035;
  const turnTube = reach * 0.008;
  const edge = reach * 0.009;
  const shaft = knobR * 0.3;
  const handles = MeshBuilder.createState(4096, 2048);
  const groups: HandleGroup[] = [];
  const group = (handle: Handle, label: string, color: Color) => {
    handles.currentGroup = groups.length;
    groups.push({ handle, label, color });
  };
  let facesCamera = false;
  switch (object.type) {
    case 'plane': {
      // the knob and the arrow slide it (moving it within itself changes nothing); the arrow points to the side it
      // cuts away
      group({ kind: 'slide', axis: LocalUp }, 'Drag to slide it along its normal', HandleColors.slide);
      knob(handles, Vec3(), knobR);
      arrow(handles, Vec3(), Vec3.create(0, (object.invert ? -1 : 1) * reach * 0.35, 0), shaft);
      group({ kind: 'turn' }, 'Drag to turn it', HandleColors.turn);
      ring(handles, reach, turnTube);
      break;
    }
    case 'sphere': {
      group({ kind: 'move' }, 'Drag to move it', HandleColors.move);
      knob(handles, Vec3(), knobR);
      group({ kind: 'scale' }, 'Drag to size it', HandleColors.outline);
      ring(handles, object.scale[0] / 2, edge * 1.4);
      facesCamera = true;
      break;
    }
    case 'cube': {
      const [sx, sy, sz] = object.scale;
      const half = [sx / 2, sy / 2, sz / 2];
      group({ kind: 'move' }, 'Drag to move it', HandleColors.move);
      knob(handles, Vec3(), knobR);
      // its faces: a knob at the center of each moves that face along its normal (the opposite one stays)
      for (const a of [0, 1, 2] as const) {
        for (const sign of [-1, 1] as const) {
          const at = Vec3();
          at[a] = sign * half[a];
          group({ kind: 'faces', faces: [{ axis: a, sign }] }, 'Drag to move this side', HandleColors.outline);
          knob(handles, at, knobR * 0.8);
        }
      }
      // the slide arrow from past the top face's knob
      group({ kind: 'slide', axis: LocalUp }, 'Drag to slide it along its up', HandleColors.slide);
      arrow(
        handles,
        Vec3.create(0, sy / 2 + knobR * 1.2, 0),
        Vec3.create(0, sy / 2 + knobR * 1.2 + reach * 0.25, 0),
        shaft,
      );
      group({ kind: 'turn' }, 'Drag to turn it', HandleColors.turn);
      ring(handles, (Math.max(sx, sz) / 2) * 1.25, turnTube);
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
            tube(handles, from, to, edge);
          }
        }
      }
      // its corners scale it whole
      group({ kind: 'scale' }, 'Drag to size it whole', HandleColors.outline);
      for (const x of [-1, 1]) {
        for (const y of [-1, 1]) {
          for (const z of [-1, 1]) knob(handles, Vec3.create(x * half[0], y * half[1], z * half[2]), edge * 2.6);
        }
      }
      break;
    }
    case 'cylinder': {
      const r = object.scale[0] / 2;
      const h = object.scale[1];
      group({ kind: 'move' }, 'Drag to move it', HandleColors.move);
      knob(handles, Vec3(), knobR);
      group({ kind: 'slide', axis: LocalUp }, 'Drag to slide it along its axis', HandleColors.slide);
      arrow(handles, Vec3.create(0, h / 2, 0), Vec3.create(0, h / 2 + reach * 0.25, 0), shaft);
      group({ kind: 'turn' }, 'Drag to turn it', HandleColors.turn);
      ring(handles, r * 1.35, turnTube);
      // its rims size its radius (about their own centers); the bottom cap its height (the top stays)
      for (const y of [-h / 2, h / 2]) {
        group(
          { kind: 'scale', axes: [0, 2], pivot: Vec3.create(0, y, 0) },
          'Drag to size its radius',
          HandleColors.outline,
        );
        ring(handles, r, edge, y);
      }
      group({ kind: 'faces', faces: [{ axis: 1, sign: -1 }] }, 'Drag to size its height', HandleColors.outline);
      knob(handles, Vec3.create(0, -h / 2, 0), edge * 2.6);
      break;
    }
    case 'infiniteCone': {
      // drawn over what it cuts: from its apex as far as their reach, along its opening (local −Z)
      const length = reach;
      const rim = length * Math.tan((coneAngle(object.scale) * Math.PI) / 180);
      group({ kind: 'move' }, 'Drag to move its apex', HandleColors.move);
      knob(handles, Vec3(), knobR);
      group({ kind: 'slide', axis: ConeOpening }, 'Drag to slide it along its opening', HandleColors.slide);
      arrow(handles, Vec3(), Vec3.create(0, 0, -reach * 0.3), shaft);
      group({ kind: 'turn' }, 'Drag to turn it', HandleColors.turn);
      ring(handles, reach * 0.12, turnTube, 0, Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, -reach * 0.15)));
      // its base rim widens or narrows it (about the base's center)
      group(
        { kind: 'custom', id: ConeAngleHandle, pivot: Vec3.create(0, 0, -length) },
        'Drag to widen or narrow it',
        HandleColors.outline,
      );
      ring(handles, rim, edge, 0, Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, -length)));
      break;
    }
  }
  return {
    // the walls as the clip volumes shape draws them
    body: clipWallsMesh(object, reach),
    handles: MeshBuilder.getMesh(handles),
    groups,
    transform: clipObjectFrame(object, bounds),
    facesCamera,
  };
}
