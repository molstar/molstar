/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import type { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Clip } from '@molstar/core/util/clip';
import { addCylinder } from './geometry/mesh/builder/cylinder.js';
import { addSphere } from './geometry/mesh/builder/sphere.js';
import type { Mesh } from './geometry/mesh/mesh.js';
import { MeshBuilder } from './geometry/mesh/mesh-builder.js';
import { Box } from './primitive/box.js';
import { Circle } from './primitive/circle.js';
import { Torus } from './primitive/torus.js';

// The volumes clip objects cut away, as meshes: their walls (translucent where drawn: two-faced, so both sides draw as
// fronts) and their edges (tubes; a box's corners rounded). Built in the object's own frame (`clipObjectFrame` places
// it in the world): a plane is a disk and a cone a cone of the given reach, as far as what they cut, being infinite.
// Shared by the clipping gizmos' bodies and the clip volumes shape (plugin `ClipVolumesShape3D`).

export type ClipObject = Clip.Props['objects'][number];

/** A clip object's local up: a plane's normal, a box's or a cylinder's axis; a cone opens along −Z */
const Up = Vec3.create(0, 1, 0);

/**
 * The object's frame in the world: its rotation, from its position; a plane's from where it crosses `around` (what it
 * cuts), so its disk is drawn over that, not at its point (which may be anywhere on the plane). The shader tests a
 * point as its `transform` maps it (MolViewSpec's clips in a grid's or a cell's space): the frame is that space's, taken
 * back to the world by the transform's inverse.
 */
export function clipObjectFrame(object: ClipObject, around: Sphere3D): Mat4 {
  const rotation = Mat4.fromRotation(Mat4(), (object.rotation.angle * Math.PI) / 180, object.rotation.axis);
  const transform = object.transform && !Mat4.isIdentity(object.transform) ? object.transform : undefined;
  let origin = Vec3.clone(object.position);
  if (object.type === 'plane') {
    const n = Vec3.transformDirection(Vec3(), Up, rotation);
    Vec3.normalize(n, n);
    const center = transform ? Vec3.transformMat4(Vec3(), around.center, transform) : around.center;
    origin = Vec3.scaleAndAdd(Vec3(), center, n, -Vec3.dot(n, Vec3.sub(Vec3(), center, object.position)));
  }
  const frame = Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), origin), rotation);
  return transform ? Mat4.mul(frame, Mat4.invert(Mat4(), transform), frame) : frame;
}

/** A cone's half-angle (degrees) from its scale (the shader's `atan2(scale.y, scale.x)`) */
export const clipConeAngle = (scale: Vec3) => (Math.atan2(scale[1], scale[0]) * 180) / Math.PI;

/** The walls of the volume the object cuts away, in its frame: a disk, a sphere, a box, a capped cylinder, a cone */
export function clipWallsMesh(object: ClipObject, reach: number, detail = 3): Mesh {
  const state = MeshBuilder.createState(1024, 512);
  state.currentGroup = 0;
  switch (object.type) {
    case 'plane':
      // `Circle` lies in XZ (its normal +Y, the plane's)
      MeshBuilder.addPrimitive(state, Mat4.fromScaling(Mat4(), Vec3.create(reach, 1, reach)), Circle({ segments: 96 }));
      break;
    case 'sphere':
      addSphere(state, Vec3(), object.scale[0] / 2, detail);
      break;
    case 'cube':
      MeshBuilder.addPrimitive(state, Mat4.fromScaling(Mat4(), Vec3.clone(object.scale)), Box());
      break;
    case 'cylinder': {
      const r = object.scale[0] / 2;
      const h = object.scale[1];
      addCylinder(state, Vec3.create(0, -h / 2, 0), Vec3.create(0, h / 2, 0), 1, {
        radiusTop: r,
        radiusBottom: r,
        radialSegments: 48,
        topCap: true,
        bottomCap: true,
      });
      break;
    }
    case 'infiniteCone': {
      const rim = reach * Math.tan((clipConeAngle(object.scale) * Math.PI) / 180);
      addCylinder(state, Vec3.create(0, 0, -reach), Vec3(), 1, { radiusTop: rim, radiusBottom: 0, radialSegments: 48 });
      // its base filled: `Circle` lies in XZ, turned to face along Z
      const base = Mat4.mul(
        Mat4(),
        Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, -reach)), Mat4.rotX90 as Mat4),
        Mat4.fromScaling(Mat4(), Vec3.create(rim, 1, rim)),
      );
      MeshBuilder.addPrimitive(state, base, Circle({ radius: 1, segments: 96 }));
      break;
    }
  }
  return twoFaced(MeshBuilder.getMesh(state));
}

/**
 * The edges of that volume, in its frame, as tubes of radius `width`: a plane's rim, three great circles of a sphere,
 * a box's twelve edges with round corners, a cylinder's rims, a cone's base rim and four lines to its apex
 */
export function clipEdgesMesh(object: ClipObject, reach: number, width: number): Mesh {
  const state = MeshBuilder.createState(2048, 1024);
  state.currentGroup = 0;
  const ring = (radius: number, t: Mat4) =>
    MeshBuilder.addPrimitive(state, t, Torus({ radius, tube: width, radialSegments: 8, tubularSegments: 64 }));
  // a ring in the frame's XZ plane (around its up), at height `y`
  const ringXZ = (radius: number, y = 0) =>
    ring(radius, Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), Vec3.create(0, y, 0)), Mat4.rotX90 as Mat4));
  const tube = (a: Vec3, b: Vec3) =>
    addCylinder(state, a, b, 1, { radiusTop: width, radiusBottom: width, radialSegments: 12 });
  switch (object.type) {
    case 'plane':
      ringXZ(reach);
      break;
    case 'sphere': {
      const r = object.scale[0] / 2;
      ringXZ(r);
      ring(r, Mat4.identity());
      ring(r, Mat4.rotY90 as Mat4);
      break;
    }
    case 'cube': {
      const half = [object.scale[0] / 2, object.scale[1] / 2, object.scale[2] / 2];
      for (const a of [0, 1, 2]) {
        const [b, c] = [0, 1, 2].filter((x) => x !== a);
        for (const sb of [-1, 1]) {
          for (const sc of [-1, 1]) {
            const from = Vec3();
            const to = Vec3();
            from[a] = -half[a];
            to[a] = half[a];
            from[b] = to[b] = sb * half[b];
            from[c] = to[c] = sc * half[c];
            tube(from, to);
          }
        }
      }
      for (const x of [-1, 1]) {
        for (const y of [-1, 1]) {
          for (const z of [-1, 1]) addSphere(state, Vec3.create(x * half[0], y * half[1], z * half[2]), width * 1.6, 2);
        }
      }
      break;
    }
    case 'cylinder': {
      const r = object.scale[0] / 2;
      const h = object.scale[1];
      ringXZ(r, -h / 2);
      ringXZ(r, h / 2);
      break;
    }
    case 'infiniteCone': {
      const rim = reach * Math.tan((clipConeAngle(object.scale) * Math.PI) / 180);
      const base = Vec3.create(0, 0, -reach);
      ring(rim, Mat4.fromTranslation(Mat4(), base));
      for (const [x, y] of [
        [1, 0],
        [-1, 0],
        [0, 1],
        [0, -1],
      ]) {
        tube(Vec3(), Vec3.create(x * rim, y * rim, -reach));
      }
      break;
    }
  }
  return MeshBuilder.getMesh(state);
}

/**
 * Each triangle twice, once facing each way: both sides of a translucent wall draw as its front. A double-sided mesh
 * would draw its back as the inside, darkened (the renderer's interior shading), which on a dark canvas is not seen.
 * Meant for walls that ignore light: the normals (one a triangle) do not show.
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
