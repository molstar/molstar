/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author ReliaSolve <russ@reliasolve.com>
 */

import { RuntimeContext, Task } from '../../mol-task';
import { ShapeProvider } from '../../mol-model/shape/provider';
import { Color } from '../../mol-util/color';
import { Kinemage } from '../../mol-io/reader/kin/schema';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Shape } from '../../mol-model/shape';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';

/// @todo Fill in geometry and coloring information

export type KinData = {
    source: Kinemage,
    transforms?: Mat4[],
}

function createKinShapeParams(kinemage?: Kinemage) {

    return {
        ...Mesh.Params,
    };
}

export const KinShapeParams = createKinShapeParams();
export type KinShapeParams = typeof KinShapeParams

function makeShapeGetter() {

    const getShape = async (ctx: RuntimeContext, kinData: KinData, props: PD.Values<KinShapeParams>, shape?: Shape<Mesh>) => {
        console.log(`XXX Number of vector lists: ${kinData.source.vectorLists.length}, ballLists: ${kinData.source.ballLists.length}, ribbonLists: ${kinData.source.ribbonLists.length}`);
        /// @todo
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
    };
    return getShape;
}

export function shapeFromKin(source: Kinemage, params?: { transforms?: Mat4[] }) {
    return Task.create<ShapeProvider<KinData, Mesh, KinShapeParams>>('Shape Provider', async ctx => {
        return {
            label: 'Mesh',
            data: { source, transforms: params?.transforms },
            params: createKinShapeParams(source),
            getShape: makeShapeGetter(),
            geometryUtils: Mesh.Utils
        };
    });
}