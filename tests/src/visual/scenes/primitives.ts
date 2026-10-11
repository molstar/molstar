/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Cylinders } from '@molstar/graphics/geo/geometry/cylinders/cylinders';
import { CylindersBuilder } from '@molstar/graphics/geo/geometry/cylinders/cylinders-builder';
import { Geometry } from '@molstar/graphics/geo/geometry/geometry';
import { Lines } from '@molstar/graphics/geo/geometry/lines/lines';
import { LinesBuilder } from '@molstar/graphics/geo/geometry/lines/lines-builder';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { Points } from '@molstar/graphics/geo/geometry/points/points';
import { PointsBuilder } from '@molstar/graphics/geo/geometry/points/points-builder';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { SpheresBuilder } from '@molstar/graphics/geo/geometry/spheres/spheres-builder';
import { Text } from '@molstar/graphics/geo/geometry/text/text';
import { TextBuilder } from '@molstar/graphics/geo/geometry/text/text-builder';
import { Box, BoxCage } from '@molstar/graphics/geo/primitive/box';
import { DodecahedronCage } from '@molstar/graphics/geo/primitive/dodecahedron';
import { Sphere } from '@molstar/graphics/geo/primitive/sphere';
import { SpikedBall } from '@molstar/graphics/geo/primitive/spiked-ball';
import { Torus } from '@molstar/graphics/geo/primitive/torus';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { packIntToRGBArray } from '@molstar/core/util/number-packing';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { TextureMesh } from '@molstar/graphics/geo/geometry/texture-mesh/texture-mesh';
import * as Shape from '@molstar/graphics/geo/shape/shape';
import { Representation } from '@molstar/graphics/repr/representation';
import { ShapeRepresentation } from '@molstar/graphics/repr/shape/representation';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

export const Palette = [
  Color(0xe41a1c),
  Color(0x377eb8),
  Color(0x4daf4a),
  Color(0x984ea3),
  Color(0xff7f00),
  Color(0xffd92f),
  Color(0xa65628),
  Color(0xf781bf),
];

export interface ShapeOptions<G extends Geometry> {
  colors?: readonly Color[];
  size?: number;
  transforms?: Mat4[];
  props?: Partial<PD.Values<Geometry.Params<G>>>;
}

/** Per-group colors (cycling through `colors`) and optional instancing via `transforms` */
export async function shapeRepr<G extends Geometry>(
  name: string,
  geometry: G,
  options: ShapeOptions<G> = {},
): Promise<Representation.Any> {
  const colors = options.colors ?? Palette;
  const size = options.size ?? 1;
  const repr = ShapeRepresentation<{}, G, Geometry.Params<G>>(
    () =>
      Shape.create(
        name,
        {},
        geometry,
        (g) => colors[g % colors.length],
        () => size,
        () => name,
        options.transforms,
      ),
    Geometry.getUtils(geometry),
  );
  await repr.createOrUpdate(options.props ?? {}, {}).run();
  return repr as unknown as Representation.Any;
}

export function rowTransforms(count: number, spacing: number): Mat4[] {
  const transforms: Mat4[] = [];
  const offset = ((count - 1) * spacing) / 2;
  for (let i = 0; i < count; ++i) {
    transforms.push(Mat4.fromTranslation(Mat4(), Vec3.create(i * spacing - offset, 0, 0)));
  }
  return transforms;
}

function translateScale(x: number, y: number, z: number, s: number) {
  const t = Mat4.fromTranslation(Mat4(), Vec3.create(x, y, z));
  return Mat4.scaleUniformly(t, t, s);
}

/** Sphere, box, torus and spiked ball in a 2x2 layout, one group each */
export function meshPrimitives(): Mesh {
  const state = MeshBuilder.createState();
  state.currentGroup = 0;
  MeshBuilder.addPrimitive(state, translateScale(-2, 2, 0, 1.6), Sphere(3));
  state.currentGroup = 1;
  MeshBuilder.addPrimitive(state, translateScale(2, 2, 0, 1.6), Box());
  state.currentGroup = 2;
  MeshBuilder.addPrimitive(
    state,
    translateScale(-2, -2, 0, 1.2),
    Torus({ tubularSegments: 48, radialSegments: 24, tube: 0.35 }),
  );
  state.currentGroup = 3;
  MeshBuilder.addPrimitive(state, translateScale(2, -2, 0, 1), SpikedBall(3));
  return MeshBuilder.getMesh(state);
}

export function meshBox(size: number): Mesh {
  const state = MeshBuilder.createState();
  state.currentGroup = 0;
  MeshBuilder.addPrimitive(state, translateScale(0, 0, 0, size), Box());
  return MeshBuilder.getMesh(state);
}

/** Thin slab behind the xy-plane, e.g. to receive shadows and occlusion */
export function meshBackdrop(size: number, z: number): Mesh {
  const state = MeshBuilder.createState();
  state.currentGroup = 0;
  const t = Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, z));
  Mat4.scale(t, t, Vec3.create(size, size, 0.2));
  MeshBuilder.addPrimitive(state, t, Box());
  return MeshBuilder.getMesh(state);
}

export function meshTorus(radius: number): Mesh {
  const state = MeshBuilder.createState();
  state.currentGroup = 0;
  MeshBuilder.addPrimitive(
    state,
    translateScale(0, 0, 0, radius),
    Torus({ tubularSegments: 64, radialSegments: 32, tube: 0.3 }),
  );
  return MeshBuilder.getMesh(state);
}

/** n x n x n grid of spheres, one group per sphere */
export function sphereGrid(n: number, spacing: number): Spheres {
  const builder = SpheresBuilder.create(n * n * n, 1);
  const offset = ((n - 1) * spacing) / 2;
  let g = 0;
  for (let i = 0; i < n; ++i) {
    for (let j = 0; j < n; ++j) {
      for (let k = 0; k < n; ++k) {
        builder.add(i * spacing - offset, j * spacing - offset, k * spacing - offset, g++);
      }
    }
  }
  return builder.getSpheres();
}

/** Cylinders from the origin along 8 directions with alternating caps */
export function cylinderStar(length: number): Cylinders {
  const builder = CylindersBuilder.create(8, 1);
  for (let i = 0; i < 8; ++i) {
    const a = (i * Math.PI) / 4;
    const x = Math.cos(a) * length,
      y = Math.sin(a) * length,
      z = (i % 2 === 0 ? 1 : -1) * length * 0.3;
    builder.add(0, 0, 0, x, y, z, 1, i % 2 === 0, true, 2, i);
  }
  return builder.getCylinders();
}

/** Points on a Fibonacci sphere */
export function pointSphere(count: number, radius: number): Points {
  const builder = PointsBuilder.create(count, 1);
  const golden = Math.PI * (3 - Math.sqrt(5));
  for (let i = 0; i < count; ++i) {
    const y = 1 - (i / (count - 1)) * 2;
    const r = Math.sqrt(1 - y * y);
    const t = golden * i;
    builder.add(Math.cos(t) * r * radius, y * radius, Math.sin(t) * r * radius, i % 8);
  }
  return builder.getPoints();
}

export function lineCages(): Lines {
  const builder = LinesBuilder.create();
  builder.addCage(translateScale(-2, 0, 0, 1.5), DodecahedronCage(), 0);
  builder.addCage(translateScale(2, 0, 0, 2.5), BoxCage(), 1);
  builder.addFixedCountDashes(Vec3.create(-4, -3, 0), Vec3.create(4, -3, 0), 12, 2);
  return builder.getLines();
}

export function textLabels(props: Partial<PD.Values<Text.Params>>): Text {
  const builder = TextBuilder.create(props, 4, 1);
  builder.add('Mol*', 0, 2, 0, 1, 1, 0);
  builder.add('Visual', -2, 0, 0, 1, 1, 1);
  builder.add('Test', 2, 0, 0, 1, 1, 2);
  builder.add('Ã…Ã„Ã– 123', 0, -2, 0, 1, 1, 3);
  return builder.getText();
}

/**
 * The triangles of `mesh` as a texture mesh (unindexed, positions, normals and packed groups in float textures),
 * as GPU-computed surfaces are, built on the CPU
 */
export function textureMeshFromMesh(graphics: GraphicsContext, mesh: Mesh): TextureMesh {
  const positions = mesh.vertexBuffer.ref.value;
  const normals = mesh.normalBuffer.ref.value;
  const groups = mesh.groupBuffer.ref.value;
  const indices = mesh.indexBuffer.ref.value;
  const vertexCount = mesh.triangleCount * 3;
  const width = Math.ceil(Math.sqrt(vertexCount));
  const height = Math.ceil(vertexCount / width);
  const position = new Float32Array(width * height * 4);
  const normal = new Float32Array(width * height * 4);
  const group = new Float32Array(width * height * 4);
  const rgb = [0, 0, 0];
  for (let i = 0; i < vertexCount; ++i) {
    const v = indices[i];
    for (let c = 0; c < 3; ++c) {
      position[i * 4 + c] = positions[v * 3 + c];
      normal[i * 4 + c] = normals[v * 3 + c];
    }
    packIntToRGBArray(groups[v], rgb, 0);
    for (let c = 0; c < 3; ++c) group[i * 4 + c] = rgb[c] / 255;
  }
  const texture = (array: Float32Array) => {
    const t = graphics.resources.texture('image-float32', 'rgba', 'float', 'nearest');
    t.load({ array, width, height });
    return t;
  };
  const groupCount = Math.max(...Array.from(groups)) + 1;
  return TextureMesh.create(
    vertexCount,
    groupCount,
    texture(position),
    texture(group),
    texture(normal),
    Sphere3D.clone(mesh.boundingSphere),
  );
}
