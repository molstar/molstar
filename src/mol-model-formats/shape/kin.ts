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
    //transparentBackfaces: PD.Select('on', PD.arrayToOptions(['off', 'on', 'opaque'] as const)),
    //doubleSided: PD.Boolean(true), // make mesh double-sided by default
    //ignoreLight: PD.Boolean(true), // ignore lighting so front/back show same color
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
  const colors: Color[] = [];

  // Every dot is in its own group because they may have colors and we look that up by group.
  let index = 0;

  for (let i = 0; i < dotLists.length; i++) {
    const positionArray = dotLists[i].positionArray;
    const colorArray = dotLists[i].colorArray;

    /// @todo Update in chunks of 100000 like the Ply files do rather than all at once like we do here.

    const numDots = positionArray.length / 3
    for (let j = 0; j < numDots; j++) {
      let group = index++;
      builderState.add(positionArray[3 * j + 0], positionArray[3 * j + 1], positionArray[3 * j + 2], group);
      // colorArray may be undefined; push a default color when not provided
      colors.push(colorArray && colorArray.length > j * 3 ?
        Color.fromRgb(255 * (colorArray[3 * j + 0]), 255 * (colorArray[3 * j + 1]), 255 * (colorArray[3 * j + 2]))
        : Color.fromRgb(255, 255, 255));
    }
  }

  const points = builderState.getPoints();
  return { points, colors };
}

async function getLines(ctx: RuntimeContext, vectorLists: VectorList[]) {
  const builderState = LinesBuilder.create();
  const widths: number[] = [];
  const colors: Color[] = [];

  // Every line is in its own group because they may have individual widths and we look
  // up the width based on the group is in the size function.
  let index = 0;

  for (let i = 0; i < vectorLists.length; i++) {
    const vertices = vectorLists[i];
    const position1Array = vertices.position1Array;
    const position2Array = vertices.position2Array;
    const widthArray = vertices.width;
    const color1Array = vertices.color1Array;
    const color2Array = vertices.color2Array;

    /// @todo Update in chunks of 100000 like the Ply files do rather than all at once like we do here.

    const numLines = position1Array.length / 3
    for (let j = 0; j < numLines; j++) {
      let group = index++;
      builderState.add(position1Array[3 * j + 0], position1Array[3 * j + 1], position1Array[3 * j + 2],
        position2Array[3 * j + 0], position2Array[3 * j + 1], position2Array[3 * j + 2],
        group);
      // widthArray may be undefined; push NaN when width not provided
      widths.push(widthArray && widthArray.length > j ? widthArray[j] : NaN);
      // colorArray may be undefined; push a default color when not provided
      colors.push(color1Array && color1Array.length > j * 3 ?
        Color.fromRgb(255 * (color1Array[3 * j + 0] + color2Array[3 * j + 0]) / 2,
                      255 * (color1Array[3 * j + 1] + color2Array[3 * j + 1]) / 2,
                      255 * (color1Array[3 * j + 2] + color2Array[3 * j + 2]) / 2)
        : Color.fromRgb(255, 255, 255));
    }
  }

  const lines = builderState.getLines();
  return { lines, widths: new Float32Array(widths), colors };
}

function addOffsetTriangle(builderState: MeshBuilder.State, a: Vec3, b: Vec3, c: Vec3, n: Vec3, offset: number) {
  const aOffset = Vec3.add(Vec3(), a, Vec3.scale(Vec3(), n, offset));
  const bOffset = Vec3.add(Vec3(), b, Vec3.scale(Vec3(), n, offset));
  const cOffset = Vec3.add(Vec3(), c, Vec3.scale(Vec3(), n, offset));
  MeshBuilder.addTriangleWithNormal(builderState, aOffset, bOffset, cOffset, n);
}

async function getMesh(ctx: RuntimeContext, ribbonObjects: RibbonObject[]) {
  const builderState = MeshBuilder.createState();
  const colors: Color[] = [];

  // Every triangle is in its own group because they may have individual colors and we look
  // up the color based on the group is in the color function.
  let index = 0;

  for (let ri = 0; ri < ribbonObjects.length; ri++) {
    const ribbonObject = ribbonObjects[ri];
    builderState.currentGroup = ri;  /// @todo Base this on something in the file instead?

    for (let i = 0; i < ribbonObject.masterArray.length; i++) {
      const coords = ribbonObject.positionArray;
      const colorArray = ribbonObject.colorArray;

      /// @todo Update in chunks of 100000 like the Ply files do rather than all at once like we do here.

      // The positionArray contains 3x as many entries as there are vertices since it's a catenation of x, y, z for each vertex.
      // There are three vertices per triangle.
      /// @todo Ribbon lighting is to be set up to make each pair of triangles look like a quad with the same normal.
      const numTriangles = coords.length / 9;
      for (let i = 0; i < numTriangles; i++) {
        const vertexList: Vec3[] = [];

        // Get the vertices for the triangle out of the position array and push them onto a list.
        for (let j = 0; j < 3; j++) {
          const v = Vec3.zero();
          v[0] = coords[3 * (3 * i + j) + 0];
          v[1] = coords[3 * (3 * i + j) + 1];
          v[2] = coords[3 * (3 * i + j) + 2];
          vertexList.push(v);
        }

        // Set the group per triangle so that we can do per-triangle coloring.
        let group = index++;
        builderState.currentGroup = group;

        // colorArray may be undefined; push a default color when not provided.
        // There is one color per group, even if we have two triangles in this group.
        /// @todo Consider averaging the colors the vertices because we can't color different vertices differently.
        const color = colorArray && colorArray.length > i * 9 ?
          Color.fromRgb(255 * colorArray[9 * i + 0],
                        255 * colorArray[9 * i + 1],
                        255 * colorArray[9 * i + 2])
          : Color.fromRgb(255, 255, 255);
        colors.push(color);

        // Find the vertics and normal for the triangle.
        let a: Vec3 = vertexList[0];
        let b: Vec3 = vertexList[1];
        let c: Vec3 = vertexList[2];
        // Flip the winding for every other triangle to keep the faces consistent.
        if (i % 2 === 1) {
          const temp = b;
          b = c;
          c = temp;
        }
        const n = Vec3.zero();

        // Put both orientations of the triangle. Add a small amount along the normal to make them
        // not be exactly on top of each other so that we only see the front face of each.
        Vec3.triangleNormal(n, a, b, c);
        addOffsetTriangle(builderState, a, b, c, n, 0.01);

        Vec3.triangleNormal(n, a, c, b);
        addOffsetTriangle(builderState, a, c, b, n, 0.01);
      }
    }
  }

  const mesh = MeshBuilder.getMesh(builderState);
  return { mesh, colors };
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
    // Get our points, adding them from all of the entries in the dot lists
    const { points: _points, colors } = await getPoints(ctx, kinData.source.dotLists);

    // Color function signature: (groupId: number, instanceId: number) => Color
    // For Lines the groupId corresponds to the line index (order added).
    const colorFn = (group: number, instance: number) => {
      return colors[group];
    }

    let _shape: Shape<Points>;
    _shape = Shape.create<Points>(
      'kin-points',
      kinData.source,
      _points,
      colorFn,                // color function reads per-point colors
      () => 1,                // size function
      () => ''                // @todo label function
    );
    return _shape;
  };
  return getShape;
}

function makeLineShapeGetter() {

    const getShape = async (ctx: RuntimeContext, kinData: KinData, props: PD.Values<KinShapeLinesParams>, shape?: Shape<Lines>) => {
        // Get our lines, adding them from all of the entries in the vector lists
        const { lines: _lines, widths, colors } = await getLines(ctx, kinData.source.vectorLists);

        // Size function signature: (groupId: number, instanceId: number) => number
        // For Lines the groupId corresponds to the line index (order added).
        const sizeFn = (group: number, instance: number) => {
          // We're specifying the radius, which is half the width.
          let w = widths[group] / 2.0;
          if (w < 1.0) { w = 1.0; }
          return Number.isFinite(w) ? w : 1.0;
        }

        // Color function signature: (groupId: number, instanceId: number) => Color
        // For Lines the groupId corresponds to the line index (order added).
        const colorFn = (group: number, instance: number) => {
          return colors[group];
        }

        let _shape: Shape<Lines>;
        _shape = Shape.create<Lines>(
          'kin-lines',
          kinData.source,
          _lines,
          colorFn,                // color function reads per-line colors
          sizeFn,                 // size function reads per-line widths
          () => ''                // @todo label function
        );
        return _shape;
    };
    return getShape;
}

function makeMeshShapeGetter() {

  const getShape = async (ctx: RuntimeContext, kinData: KinData, props: PD.Values<KinShapeMeshParams>, shape?: Shape<Mesh>) => {

    let { mesh: _mesh, colors } = await getMesh(ctx, kinData.source.ribbonLists);
    // Ensure that _mesh is not undifined before we pass it to Shape.create.  If it is undefined, create an empty mesh instead.
    if (!_mesh) {
      console.warn('No mesh could be created from the KIN data.  Creating an empty mesh instead.');
      _mesh = Mesh.createEmpty();
    }

    // Color function signature: (groupId: number, instanceId: number) => Color
    // For Lines the groupId corresponds to the line index (order added).
    const colorFn = (group: number, instance: number) => {
      return colors[group];
    }

    let _shape: Shape<Mesh>;
    _shape = Shape.create<Mesh>(
      'kin-mesh',
      kinData.source,
      _mesh,
      colorFn,                // color function reads per-triangle colors
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
