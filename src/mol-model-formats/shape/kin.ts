/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author ReliaSolve <russ@reliasolve.com>
 */

import { RuntimeContext, Task } from '../../mol-task';
import { ShapeProvider } from '../../mol-model/shape/provider';
import { Color } from '../../mol-util/color';
import { Kinemage, DotList, VectorList, RibbonObject, BallList } from '../../mol-io/reader/kin/schema';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../mol-geo/geometry/lines/lines-builder';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Points } from '../../mol-geo/geometry/points/points';
import { PointsBuilder } from '../../mol-geo/geometry/points/points-builder';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { SpheresBuilder } from '../../mol-geo/geometry/spheres/spheres-builder';
import { Shape } from '../../mol-model/shape';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
//import { ValueCell } from '../../mol-util/value-cell';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';

/// @todo Fill in geometry and coloring information

export type KinData = {
    source: Kinemage,
    transforms?: Mat4[],
}

function createKinShapePointsParams(kinemage?: Kinemage) {

  return {
    ...Points.Params,
  };
}
export const KinShapePointsParams = createKinShapePointsParams();
export type KinShapePointsParams = typeof KinShapePointsParams
function createKinShapeLinesParams(kinemage?: Kinemage) {

    return {
        ...Lines.Params,
    };
}
export const KinShapeLinesParams = createKinShapeLinesParams();
export type KinShapeLinesParams = typeof KinShapeLinesParams
function createKinShapeMeshParams(kinemage?: Kinemage) {

  return {
    ...Mesh.Params,
  };
}

export const KinShapeMeshParams = createKinShapeMeshParams();
export type KinShapeMeshParams = typeof KinShapeMeshParams

function createKinShapeSpheresParams(kinemage?: Kinemage) {

  return {
    ...Spheres.Params,
  };
}

export const KinShapeSpheresParams = createKinShapeSpheresParams();
export type KinShapeSpheresParams = typeof KinShapeSpheresParams

async function getPoints(ctx: RuntimeContext, dotLists: DotList[]) {
  const builderState = PointsBuilder.create();

  for (let i = 0; i < dotLists.length; i++) {
    const positionArray = dotLists[i].positionArray;

    /// @todo Update in chunks of 100000 like the Ply files do rather than all at once like we do here.

    const group = i;  /// @todo Base this on something in the file instead?
    const numDots = positionArray.length / 3
    for (let j = 0; j < numDots; j++) {
      builderState.add(positionArray[3 * j + 0], positionArray[3 * j + 1], positionArray[3 * j + 2], group);
    }
  }

  return builderState.getPoints();
}

async function getLines(ctx: RuntimeContext, vectorLists: VectorList[]) {
  const builderState = LinesBuilder.create();

  for (let i = 0; i < vectorLists.length; i++) {
    const vertices = vectorLists[i];
    const position1Array = vertices.position1Array;
    const position2Array = vertices.position2Array;

    /// @todo Update in chunks of 100000 like the Ply files do rather than all at once like we do here.

    const group = i;  /// @todo Base this on something in the file instead?
    const numLines = position1Array.length / 3
    for (let j = 0; j < numLines; j++) {
      builderState.add(position1Array[3 * j + 0], position1Array[3 * j + 1], position1Array[3 * j + 2],
        position2Array[3 * j + 0], position2Array[3 * j + 1], position2Array[3 * j + 2],
        group);
    }
  }

  return builderState.getLines();
}

async function getMesh(ctx: RuntimeContext, ribbonObjects: RibbonObject[]) {
  const builderState = MeshBuilder.createState();

  for (let ri = 0; ri < ribbonObjects.length; ri++) {
    const ribbonObject = ribbonObjects[ri];
    builderState.currentGroup = ri;  /// @todo Base this on something in the file instead?

    for (let i = 0; i < ribbonObject.masterArray.length; i++) {
      const coords = ribbonObject.positionArray;

      /// @todo Update in chunks of 100000 like the Ply files do rather than all at once like we do here.

      // The positionArray contains 3x as many entries as there are vertices since it's a catenation of x, y, z for each vertex.
      // The breakArray contains a boolean for each vertex indicating if there is a break there.
      // We need to add a triangle for each new verticex after the first two, but we also need to check the break array to see
      // if there is a break between any of the vertices in the triangle.  If there is a break, we need to start over accumulating
      // vertices until we get to three.
      // Keep track of the parity and flip over every other triangle so their front faces match.
      /// @todo Lighting is to be set up to make each pair of triangles look like a quad with the same normal.
      const numVertices = coords.length / 3;
      const vertexList: Vec3[] = [];
      let flip = true;  // Start with true so that the third flip will make it false.
      for (let i = 0; i < numVertices; i++) {
        // Get the next vertex
        const v = Vec3.zero();
        v[0] = coords[3 * i + 0];
        v[1] = coords[3 * i + 1];
        v[2] = coords[3 * i + 2];
        vertexList.push(v);
        flip = !flip;

        // Once we have three vertices, make a triangle out of them and pop the oldest off the list.
        if (vertexList.length == 3) {
          if (flip) {
            MeshBuilder.addTriangle(builderState, vertexList[2], vertexList[1], vertexList[0]);
          } else {
            MeshBuilder.addTriangle(builderState, vertexList[0], vertexList[1], vertexList[2]);
          }
          vertexList.shift();
        }

        // If there is a break, clear the vertex list so we start accumulating vertices for the next triangle from scratch.
        if (ribbonObject.breakArray[i]) {
          vertexList.length = 0;
          flip = true;
        }
      }
    }
  }

  return MeshBuilder.getMesh(builderState);
}

/**
 * Build spheres geometry and collect per-sphere radii from the KIN BallList entries.
 * Returns an object with the Spheres geometry and a Float32Array with per-center radii (one entry per center, in the same order they were added).
 */
async function getSpheres(ctx: RuntimeContext, balls: BallList[]) {
  const builderState = SpheresBuilder.create();
  const radii: number[] = [];

  // Every ball is in its own group because they may have individual radii and we look
  // up the radius based on the group is in the size function.
  let index = 0;

  for (let i = 0; i < balls.length; i++) {
    const positionArray = balls[i].positionArray;
    const radiusArray = balls[i].radiusArray;
    /// @todo Update in chunks of 100000 like the Ply files do rather than all at once like we do here.

    const numBalls = positionArray.length / 3;
    for (let j = 0; j < numBalls; j++) {
      const group = index++;
      builderState.add(positionArray[3 * j + 0], positionArray[3 * j + 1], positionArray[3 * j + 2], group);
      // radiusArray may be undefined; push NaN when radius not provided
      radii.push(radiusArray && radiusArray.length > j ? radiusArray[j] : NaN);
    }
  }

  const spheres = builderState.getSpheres();
  return { spheres, radii: new Float32Array(radii) };
}

function makePointsShapeGetter() {

  const getShape = async (ctx: RuntimeContext, kinData: KinData, props: PD.Values<KinShapePointsParams>, shape?: Shape<Points>) => {
    console.log(`XXX Number of dot lists for points: ${kinData.source.dotLists.length}`);
    // Get our points, adding them from all of the entries in the dot lists
    const _points = await getPoints(ctx, kinData.source.dotLists);

    let _shape: Shape<Points>;
    _shape = Shape.create<Points>(
      'kin-points',
      kinData.source,
      _points,
      () => Color(0x7F7F7F),  // @todo color function
      () => 1,                // size function
      () => ''                // @todo label function
    );
    return _shape;
  };
  return getShape;
}

function makeLineShapeGetter() {

    const getShape = async (ctx: RuntimeContext, kinData: KinData, props: PD.Values<KinShapeLinesParams>, shape?: Shape<Lines>) => {
        console.log(`XXX Number of vector lists for lines: ${kinData.source.vectorLists.length}`);
        // Get our lines, adding them from all of the entries in the vector lists
        const _lines = await getLines(ctx, kinData.source.vectorLists);

        let _shape: Shape<Lines>;
        _shape = Shape.create<Lines>(
          'kin-lines',
          kinData.source,
          _lines,
          () => Color(0x7F7F7F),  // @todo color function
          () => 1,                // size function
          () => ''                // @todo label function
        );
        return _shape;
    };
    return getShape;
}

function makeMeshShapeGetter() {

  const getShape = async (ctx: RuntimeContext, kinData: KinData, props: PD.Values<KinShapeMeshParams>, shape?: Shape<Mesh>) => {
    console.log(`XXX Number of ribbon lists for mesh: ribbonLists: ${kinData.source.ribbonLists.length}`);
    /// @todo

    let _mesh = await getMesh(ctx, kinData.source.ribbonLists);
    // Ensure that _mesh is not undifined before we pass it to Shape.create.  If it is undefined, create an empty mesh instead.
    if (!_mesh) {
      console.warn('No mesh could be created from the KIN data.  Creating an empty mesh instead.');
      _mesh = Mesh.createEmpty();
    }

    let _shape: Shape<Mesh>;
    _shape = Shape.create<Mesh>(
      'kin-mesh',
      kinData.source,
      _mesh,
      () => Color(0x7F7F7F),  // @todo color function
      () => 1,                // size function
      () => ''                // @todo label function
    );
    return _shape;
  };
  return getShape;
}

/**
 * Spheres shape getter: uses per-center radii read from the KIN BallList radiusArray when available.
 */
function makeSpheresShapeGetter() {

  const getShape = async (ctx: RuntimeContext, kinData: KinData, props: PD.Values<KinShapeSpheresParams>, shape?: Shape<Spheres>) => {
    console.log(`XXX Number of ball lists for spheres: ${kinData.source.ballLists.length}`);
    // Build spheres geometry and collect per-center radii
    const { spheres: _spheres, radii } = await getSpheres(ctx, kinData.source.ballLists);

    // size function signature: (groupId: number, instanceId: number) => number
    // For Spheres the groupId corresponds to the center index (order added).
    const sizeFn = (group: number, instance: number) => {
      const r = radii[group];
      return Number.isFinite(r) ? r : 1.0;
    };

    let _shape: Shape<Spheres>;
    _shape = Shape.create<Spheres>(
      'kin-spheres',
      kinData.source,
      _spheres,
      () => Color(0x7F7F7F),  // @todo color function
      sizeFn,                 // size function reads per-center radii
      () => ''                // @todo label function
    );
    return _shape;
  };
  return getShape;
}

export function shapePointsFromKin(source: Kinemage, params?: { transforms?: Mat4[] }) {
  return Task.create<ShapeProvider<KinData, Points, KinShapePointsParams>>('Kin Shape Points Provider', async ctx => {
    return {
      label: 'Points',
      data: { source, transforms: params?.transforms },
      params: createKinShapePointsParams(source),
      getShape: makePointsShapeGetter(),
      geometryUtils: Points.Utils
    };
  });
}

export function shapeLinesFromKin(source: Kinemage, params?: { transforms?: Mat4[] }) {
  return Task.create<ShapeProvider<KinData, Lines, KinShapeLinesParams>>('Kin Shape Lines Provider', async ctx => {
    return {
      label: 'Lines',
      data: { source, transforms: params?.transforms },
      params: createKinShapeLinesParams(source),
      getShape: makeLineShapeGetter(),
      geometryUtils: Lines.Utils
    };
  });
}

export function shapeMeshFromKin(source: Kinemage, params?: { transforms?: Mat4[] }) {
  return Task.create<ShapeProvider<KinData, Mesh, KinShapeMeshParams>>('Kin Shape Mesh Provider', async ctx => {
    return {
      label: 'Mesh',
      data: { source, transforms: params?.transforms },
      params: createKinShapeMeshParams(source),
      getShape: makeMeshShapeGetter(),
      geometryUtils: Mesh.Utils
    };
  });
}

export function shapeSpheresFromKin(source: Kinemage, params?: { transforms?: Mat4[] }) {
  return Task.create<ShapeProvider<KinData, Spheres, KinShapeSpheresParams>>('Kin Shape Spheres Provider', async ctx => {
    return {
      label: 'Spheres',
      data: { source, transforms: params?.transforms },
      params: createKinShapeSpheresParams(source),
      getShape: makeSpheresShapeGetter(),
      geometryUtils: Spheres.Utils
    };
  });
}
