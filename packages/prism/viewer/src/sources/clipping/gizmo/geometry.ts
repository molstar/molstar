/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
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
// covers what it cuts), and its handles, opaque, one group each: what dragging it does. Built in the cut's frame (its
// rotation, from its center; a plane's from where it crosses what it cuts), so a move or a turn only moves the frame.
// Sized by what it cuts (`reach`, the radius of their bounds): knobs and rings stay in proportion to the scene.

/** The handles' groups, in order: a group's id is its index */
export const HandleGroups: readonly Handle[] = [
  'move',
  'slide',
  'turn',
  'size',
  'size-x',
  'size-y',
  'size-z',
  'radius',
];

export const HandleLabels: Readonly<Record<Handle, string>> = {
  move: 'Drag to move the cut',
  slide: 'Drag to slide it along its axis',
  turn: 'Drag to turn it',
  size: 'Drag to size it',
  'size-x': 'Drag to size its width',
  'size-y': 'Drag to size its height',
  'size-z': 'Drag to size its depth',
  radius: 'Drag to size its radius',
};

const group = (h: Handle) => HandleGroups.indexOf(h);

export interface GizmoGeometry {
  readonly body: Mesh;
  readonly handles: Mesh;
  /** The handles it has (a plane has no size, a sphere no axis or turn) */
  readonly has: readonly Handle[];
  /** The cut's frame in the world: where it sits and how it is turned */
  readonly transform: Mat4;
  /** Its axis in the world: a plane's normal, a cylinder's or a box's up, a cone's opening */
  readonly axis: Vec3;
}

const Up = Vec3.create(0, 1, 0);

/** The meshes of a cut's gizmo over what it cuts */
export function gizmoGeometry(object: CutObject, bounds: Sphere3D): GizmoGeometry {
  const reach = Math.max(bounds.radius, 1);
  const knob = reach * 0.035;
  const tube = reach * 0.008;
  const rotation = Mat4.fromRotation(Mat4(), (object.rotation.angle * Math.PI) / 180, object.rotation.axis);
  const body = MeshBuilder.createState(1024, 512);
  const handles = MeshBuilder.createState(1024, 512);
  const has: Handle[] = ['move'];
  const knobAt = (group: number, at: Vec3) => {
    handles.currentGroup = group;
    addSphere(handles, at, knob, 2);
  };
  // A ring in the frame's XZ plane (around its up), of radius `r`, at height `y`
  const ringXZ = (state: MeshBuilder.State, group: number, r: number, y = 0) => {
    state.currentGroup = group;
    const m = Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), Vec3.create(0, y, 0)), Mat4.rotX90 as Mat4);
    MeshBuilder.addPrimitive(state, m, Torus({ radius: r, tube, radialSegments: 8, tubularSegments: 64 }));
  };
  const arrow = (group: number, from: Vec3, to: Vec3) => {
    handles.currentGroup = group;
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
      const disk = Mat4.fromScaling(Mat4(), Vec3.create(reach, 1, reach));
      MeshBuilder.addPrimitive(body, disk, Circle({ radius: 1, segments: 96 }));
      // the arrow points to the side it cuts away
      arrow(1, Vec3(), Vec3.create(0, (object.invert ? -1 : 1) * reach * 0.35, 0));
      ringXZ(handles, 2, reach);
      has.push('slide', 'turn');
      break;
    }
    case 'sphere': {
      const r = object.scale[0] / 2;
      body.currentGroup = 0;
      addSphere(body, Vec3(), r, 3);
      ringXZ(handles, 3, r);
      has.push('size');
      break;
    }
    case 'cube': {
      const [sx, sy, sz] = object.scale;
      body.currentGroup = 0;
      MeshBuilder.addPrimitive(body, Mat4.fromScaling(Mat4(), Vec3.create(sx, sy, sz)), Box());
      arrow(1, Vec3.create(0, sy / 2, 0), Vec3.create(0, sy / 2 + reach * 0.25, 0));
      ringXZ(handles, 2, (Math.max(sx, sz) / 2) * 1.2);
      knobAt(group('size'), Vec3.create(sx / 2, sy / 2, sz / 2));
      knobAt(group('size-x'), Vec3.create(sx / 2, 0, 0));
      knobAt(group('size-y'), Vec3.create(0, -sy / 2, 0));
      knobAt(group('size-z'), Vec3.create(0, 0, sz / 2));
      has.push('slide', 'turn', 'size', 'size-x', 'size-y', 'size-z');
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
      arrow(1, Vec3.create(0, h / 2, 0), Vec3.create(0, h / 2 + reach * 0.25, 0));
      ringXZ(handles, 2, r * 1.25);
      knobAt(group('radius'), Vec3.create(r, 0, 0));
      knobAt(group('size-y'), Vec3.create(0, -h / 2, 0));
      has.push('slide', 'turn', 'radius', 'size-y');
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
      arrow(1, Vec3(), Vec3.create(0, 0, -reach * 0.3));
      handles.currentGroup = 2;
      MeshBuilder.addPrimitive(
        handles,
        Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, -length)),
        Torus({ radius: rim, tube, radialSegments: 8, tubularSegments: 64 }),
      );
      knobAt(3, Vec3.create(rim, 0, -length));
      has.push('slide', 'turn', 'size');
      break;
    }
  }
  // a plane's knob slides it (moving it within itself changes nothing); every other cut's moves it in the screen's plane
  knobAt(object.type === 'plane' ? 1 : 0, Vec3());
  if (object.type === 'plane') has.splice(has.indexOf('move'), 1);
  const transform = Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), origin), rotation);
  return { body: twoFaced(MeshBuilder.getMesh(body)), handles: MeshBuilder.getMesh(handles), has, transform, axis };
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
  return object.type === 'infiniteCone' ? rotated(object.rotation, ConeOpening) : rotated(object.rotation, Up);
}
