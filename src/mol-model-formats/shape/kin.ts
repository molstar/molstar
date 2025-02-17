/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author ReliaSolve <russ@reliasolve.com>
 */

import { RuntimeContext, Task } from '../../mol-task';
import { ShapeProvider } from '../../mol-model/shape/provider';
import { Color } from '../../mol-util/color';
import { Kinemage, VectorList } from '../../mol-io/reader/kin/schema';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../mol-geo/geometry/lines/lines-builder';
//import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Shape } from '../../mol-model/shape';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
//import { ValueCell } from '../../mol-util/value-cell';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';

/// @todo Fill in geometry and coloring information

export type KinData = {
    source: Kinemage,
    transforms?: Mat4[],
}

function createKinShapeParams(kinemage?: Kinemage) {

    return {
        ...Lines.Params,
    };
}

export const KinShapeParams = createKinShapeParams();
export type KinShapeParams = typeof KinShapeParams

async function getLines(ctx: RuntimeContext, vectorLists: VectorList[]) {
  const builderState = LinesBuilder.create();

  for (let i = 0; i < vectorLists.length; i++) {
    const vertices = vectorLists[i];
    const position1Array = vertices.position1Array;
    const position2Array = vertices.position2Array;

    /// @todo Update in chunks of 100000 like the Ply files do rather than all at once like we do here.

    const group = i;  /// @todo Base this on something in the file instead?
    const numLines = position1Array.length / 3
    for (let i = 0; i < numLines; i++) {
      builderState.add(position1Array[3 * i + 0], position1Array[3 * i + 1], position1Array[3 * i + 2],
        position2Array[3 * i + 0], position2Array[3 * i + 1], position2Array[3 * i + 2],
        group);

      if (ctx.shouldUpdate && (i % 10000 == 0)) {
        await ctx.update({ message: 'adding kin line vertices', current: i, max: numLines });
      }
    }
  }

  return builderState.getLines();
}

function makeShapeGetter() {

    const getShape = async (ctx: RuntimeContext, kinData: KinData, props: PD.Values<KinShapeParams>, shape?: Shape<Lines>) => {
        console.log(`XXX Number of vector lists: ${kinData.source.vectorLists.length}, ballLists: ${kinData.source.ballLists.length}, ribbonLists: ${kinData.source.ribbonLists.length}`);
        /// @todo

        /*
        // Create an empty Mesh
        const mesh = Mesh.createEmpty();

        // Create an empty Shape with the empty Mesh
        const emptyShape = Shape.create(
          'Empty Shape', // id
          kinData,      // source data
          mesh,         // geometry
          () => Color(0xFFFFFF), // color function
          () => 1,      // size function
          () => ''      // label function
        );

        return emptyShape;
        */

        // Get our lines, adding them from all of the entries in the vector lists
        const _lines = await getLines(ctx, kinData.source.vectorLists);

        /*
        let _lines: Lines = {};
        for (let i = 0; i < kinData.source.vectorLists.length; i++) {
          _lines = getLines(ctx, kinData.source.vectorLists[i], [], _lines);
        }
        */

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

export function shapeFromKin(source: Kinemage, params?: { transforms?: Mat4[] }) {
    return Task.create<ShapeProvider<KinData, Lines, KinShapeParams>>('Shape Provider', async ctx => {
        return {
            label: 'Lines',
            data: { source, transforms: params?.transforms },
            params: createKinShapeParams(source),
            getShape: makeShapeGetter(),
            geometryUtils: Lines.Utils
        };
    });
}