/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';
import { addCylinder, addSimpleCylinder } from '@molstar/graphics/geo/geometry/mesh/builder/cylinder';
import { addSphere } from '@molstar/graphics/geo/geometry/mesh/builder/sphere';
import type { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { Box } from '@molstar/graphics/geo/primitive/box';
import { Circle } from '@molstar/graphics/geo/primitive/circle';
import { Torus } from '@molstar/graphics/geo/primitive/torus';
import { ConeOpening, type CutObject, coneAngle, LocalUp, rotated } from '../cuts.js';
import type { Handle } from './gestures.js';

// A cut's gizmo (clipping.md §4): its body, translucent and not picked (a click needs alpha ≥ 0.5, and the body
// covers what it cuts), and its handles, opaque, one group each, each saying what dragging it does (`Handle`). Built
// in the cut's frame (its rotation, from its center; a plane's from where it crosses what it cuts), so a move or a
// turn only moves the frame. Sized by what it cuts (`reach`, the radius of their bounds): knobs, tubes and rings stay
// in proportion to the scene. A closed shape's outline is its handles: a box's edges and corners, a cylinder's rims,
// a cone's base rim, a sphere's ring (turned to face the camera, `facesCamera`).

/** A handle of a gizmo: what dragging it does, how it reads, its color */
export interface HandleGroup {
  readonly handle: Handle;
  readonly label: string;
  readonly color: Color;
}

/** The colors by what a handle does: the outline's (sizing) one light accent, the rest as before */
export const HandleColors = {
  move: Color(0xf5a623),
  slide: Color(0x4a90e2),
  turn: Color(0x50b860),
  outline: Color(0x9db8e0),
} as const;

export interface GizmoGeometry {
  readonly body: Mesh;
  readonly handles: Mesh;
  /** The handles' groups: a group's id is its index */
  readonly groups: readonly HandleGroup[];
  /** The cut's frame in the world: where it sits and how it is turned */
  readonly transform: Mat4;
  /** Its axis in the world: a plane's normal, a cylinder's or a box's up, a cone's opening */
  readonly axis: Vec3;
  /** Its frame turned to face the camera, not by the cut (a sphere: its ring always round, always to grab) */
  readonly facesCamera: boolean;
}

/** The meshes of a cut's gizmo over what it cuts */
export function gizmoGeometry(object: CutObject, bounds: Sphere3D): GizmoGeometry {
  const reach = Math.max(bounds.radius, 1);
  const knob = reach * 0.035;
  const tube = reach * 0.008;
  const edge = reach * 0.009;
  const rotation = Mat4.fromRotation(Mat4(), (object.rotation.angle * Math.PI) / 180, object.rotation.axis);
  const body = MeshBuilder.createState(1024, 512);
  const handles = MeshBuilder.createState(4096, 2048);
  const groups: HandleGroup[] = [];
  const group = (handle: Handle, label: string, color: Color) => {
    handles.currentGroup = groups.length;
    groups.push({ handle, label, color });
  };
  const knobAt = (at: Vec3, radius = knob) => addSphere(handles, at, radius, 2);
  const tubeBetween = (a: Vec3, b: Vec3, radius = edge) =>
    addCylinder(handles, a, b, 1, { radiusTop: radius, radiusBottom: radius, radialSegments: 12 });
  // A ring in the frame's XZ plane (around its up), of radius `r`, at height `y`
  const ringXZ = (state: MeshBuilder.State, r: number, y = 0, width = tube) => {
    const m = Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), Vec3.create(0, y, 0)), Mat4.rotX90 as Mat4);
    MeshBuilder.addPrimitive(state, m, Torus({ radius: r, tube: width, radialSegments: 8, tubularSegments: 64 }));
  };
  const arrow = (from: Vec3, to: Vec3) => {
    const shaftEnd = Vec3.lerp(Vec3(), from, to, 0.8);
    addCylinder(handles, from, shaftEnd, 1, { radiusTop: knob * 0.3, radiusBottom: knob * 0.3, radialSegments: 12 });
    // as MolViewSpec's arrows: `addSimpleCylinder` puts `radiusBottom` at its start (`addCylinder` does not)
    addSimpleCylinder(handles, shaftEnd, to, {
      radiusTop: 0,
      radiusBottom: knob * 0.9,
      radialSegments: 16,
      bottomCap: true,
    });
  };

  let origin = Vec3.clone(object.position);
  let axis = rotated(object.rotation, LocalUp);
  let facesCamera = false;
  switch (object.type) {
    case 'plane': {
      // where the plane crosses what it cuts, not its point (which keys and edits move along the plane)
      const n = axis;
      origin = Vec3.scaleAndAdd(
        Vec3(),
        bounds.center,
        n,
        -Vec3.dot(n, Vec3.sub(Vec3(), bounds.center, object.position)),
      );
      body.currentGroup = 0;
      // `Circle` lies in XZ already (its normal +Y, the plane's)
      MeshBuilder.addPrimitive(body, Mat4.fromScaling(Mat4(), Vec3.create(reach, 1, reach)), Circle({ segments: 96 }));
      // the knob and the arrow slide it (moving it within itself changes nothing); the arrow points to the side it
      // cuts away
      group({ kind: 'slide' }, 'Drag to slide it along its normal', HandleColors.slide);
      knobAt(Vec3());
      arrow(Vec3(), Vec3.create(0, (object.invert ? -1 : 1) * reach * 0.35, 0));
      group({ kind: 'turn' }, 'Drag to turn it', HandleColors.turn);
      ringXZ(handles, reach);
      break;
    }
    case 'sphere': {
      const r = object.scale[0] / 2;
      body.currentGroup = 0;
      addSphere(body, Vec3(), r, 3);
      group({ kind: 'move' }, 'Drag to move it', HandleColors.move);
      knobAt(Vec3());
      group({ kind: 'size' }, 'Drag to size it', HandleColors.outline);
      ringXZ(handles, r, 0, edge * 1.4);
      facesCamera = true;
      break;
    }
    case 'cube': {
      const [sx, sy, sz] = object.scale;
      const half = [sx / 2, sy / 2, sz / 2];
      body.currentGroup = 0;
      MeshBuilder.addPrimitive(body, Mat4.fromScaling(Mat4(), Vec3.create(sx, sy, sz)), Box());
      group({ kind: 'move' }, 'Drag to move it', HandleColors.move);
      knobAt(Vec3());
      // its faces: a knob at the center of each moves that face along its normal (the opposite one stays)
      for (const a of [0, 1, 2] as const) {
        for (const sign of [-1, 1] as const) {
          const at = Vec3();
          at[a] = sign * half[a];
          group({ kind: 'faces', faces: [{ axis: a, sign }] }, 'Drag to move this side', HandleColors.outline);
          knobAt(at, knob * 0.8);
        }
      }
      // the slide arrow from past the top face's knob
      group({ kind: 'slide' }, 'Drag to slide it along its up', HandleColors.slide);
      arrow(Vec3.create(0, sy / 2 + knob * 1.2, 0), Vec3.create(0, sy / 2 + knob * 1.2 + reach * 0.25, 0));
      group({ kind: 'turn' }, 'Drag to turn it', HandleColors.turn);
      ringXZ(handles, (Math.max(sx, sz) / 2) * 1.25);
      // its edges, the four along each axis one handle (lit together): the section across that axis scaled about the
      // center
      for (const a of [0, 1, 2] as const) {
        const [b, c] = ([0, 1, 2] as const).filter((x) => x !== a);
        group({ kind: 'across', axis: a }, 'Drag to size it across these edges', HandleColors.outline);
        for (const sb of [-1, 1]) {
          for (const sc of [-1, 1]) {
            const from = Vec3();
            const to = Vec3();
            from[a] = -half[a];
            to[a] = half[a];
            from[b] = to[b] = sb * half[b];
            from[c] = to[c] = sc * half[c];
            tubeBetween(from, to);
          }
        }
      }
      // its corners scale it whole
      group({ kind: 'size' }, 'Drag to size it whole', HandleColors.outline);
      for (const x of [-1, 1]) {
        for (const y of [-1, 1])
          for (const z of [-1, 1]) knobAt(Vec3.create(x * half[0], y * half[1], z * half[2]), edge * 2.6);
      }
      break;
    }
    case 'cylinder': {
      const r = object.scale[0] / 2;
      const h = object.scale[1];
      body.currentGroup = 0;
      addCylinder(body, Vec3.create(0, -h / 2, 0), Vec3.create(0, h / 2, 0), 1, {
        radiusTop: r,
        radiusBottom: r,
        radialSegments: 48,
        topCap: true,
        bottomCap: true,
      });
      group({ kind: 'move' }, 'Drag to move it', HandleColors.move);
      knobAt(Vec3());
      group({ kind: 'slide' }, 'Drag to slide it along its axis', HandleColors.slide);
      arrow(Vec3.create(0, h / 2, 0), Vec3.create(0, h / 2 + reach * 0.25, 0));
      group({ kind: 'turn' }, 'Drag to turn it', HandleColors.turn);
      ringXZ(handles, r * 1.35);
      // its rims size its radius (about their own centers); the bottom cap its height (the top stays)
      for (const y of [-h / 2, h / 2]) {
        group({ kind: 'radius', pivot: Vec3.create(0, y, 0) }, 'Drag to size its radius', HandleColors.outline);
        ringXZ(handles, r, y, edge);
      }
      group({ kind: 'faces', faces: [{ axis: 1, sign: -1 }] }, 'Drag to size its height', HandleColors.outline);
      knobAt(Vec3.create(0, -h / 2, 0), edge * 2.6);
      break;
    }
    case 'infiniteCone': {
      // drawn over what it cuts: from its apex as far as their reach, along its opening (local −Z)
      axis = rotated(object.rotation, ConeOpening);
      const length = reach;
      const rim = length * Math.tan((coneAngle(object.scale) * Math.PI) / 180);
      body.currentGroup = 0;
      addCylinder(body, Vec3.create(0, 0, -length), Vec3(), 1, { radiusTop: rim, radiusBottom: 0, radialSegments: 48 });
      // its base filled: `Circle` lies in XZ, turned to face along Z
      const base = Mat4.mul(
        Mat4(),
        Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, -length)), Mat4.rotX90 as Mat4),
        Mat4.fromScaling(Mat4(), Vec3.create(rim, 1, rim)),
      );
      MeshBuilder.addPrimitive(body, base, Circle({ radius: 1, segments: 96 }));
      group({ kind: 'move' }, 'Drag to move its apex', HandleColors.move);
      knobAt(Vec3());
      group({ kind: 'slide' }, 'Drag to slide it along its opening', HandleColors.slide);
      arrow(Vec3(), Vec3.create(0, 0, -reach * 0.3));
      group({ kind: 'turn' }, 'Drag to turn it', HandleColors.turn);
      MeshBuilder.addPrimitive(
        handles,
        Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, -reach * 0.15)),
        Torus({ radius: reach * 0.12, tube, radialSegments: 8, tubularSegments: 48 }),
      );
      // its base rim widens or narrows it (about the base's center)
      group({ kind: 'angle', pivot: Vec3.create(0, 0, -length) }, 'Drag to widen or narrow it', HandleColors.outline);
      MeshBuilder.addPrimitive(
        handles,
        Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, -length)),
        Torus({ radius: rim, tube: edge, radialSegments: 8, tubularSegments: 64 }),
      );
      break;
    }
  }
  const transform = Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), origin), rotation);
  return {
    body: twoFaced(MeshBuilder.getMesh(body)),
    handles: MeshBuilder.getMesh(handles),
    groups,
    transform,
    axis,
    facesCamera,
  };
}

/** A frame at `transform`'s place turned so its up (+Y) points from there to the camera (`eye`) */
export function facingCamera(transform: Mat4, eye: Vec3): Mat4 {
  const at = Vec3.fromArray(Vec3(), transform, 12);
  const toEye = Vec3.sub(Vec3(), eye, at);
  if (Vec3.magnitude(toEye) < 1e-9) return transform;
  Vec3.normalize(toEye, toEye);
  const turn = Mat4.fromRotation(Mat4(), Math.acos(Math.max(-1, Math.min(1, Vec3.dot(LocalUp, toEye)))), axisTo(toEye));
  return Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), at), turn);
}

function axisTo(to: Vec3): Vec3 {
  const a = Vec3.cross(Vec3(), LocalUp, to);
  return Vec3.magnitude(a) < 1e-6 ? Vec3.create(1, 0, 0) : Vec3.normalize(a, a);
}

/**
 * Each triangle twice, once facing each way: both sides of a translucent body draw as its front. A double-sided mesh
 * would draw its back as the inside, darkened (the renderer's interior shading), which on a dark canvas is not seen.
 * The body ignores light, so its normals (one a triangle) do not show.
 */
export function twoFaced(mesh: Mesh): Mesh {
  const v = mesh.vertexBuffer.ref.value;
  const idx = mesh.indexBuffer.ref.value;
  const groups = mesh.groupBuffer.ref.value;
  const count = mesh.triangleCount;
  const out = MeshBuilder.createState(count * 6, count * 2);
  const a = Vec3();
  const b = Vec3();
  const c = Vec3();
  for (let t = 0; t < count; t++) {
    const i = idx[3 * t];
    Vec3.fromArray(a, v, 3 * i);
    Vec3.fromArray(b, v, 3 * idx[3 * t + 1]);
    Vec3.fromArray(c, v, 3 * idx[3 * t + 2]);
    out.currentGroup = groups[i];
    MeshBuilder.addTriangle(out, a, b, c);
    MeshBuilder.addTriangle(out, a, c, b);
  }
  return MeshBuilder.getMesh(out);
}

/** The axis a slide follows for a cut: the plane's normal, the box's or cylinder's up, the cone's opening */
export function slideAxis(object: CutObject): Vec3 {
  return object.type === 'infiniteCone' ? rotated(object.rotation, ConeOpening) : rotated(object.rotation, LocalUp);
}
