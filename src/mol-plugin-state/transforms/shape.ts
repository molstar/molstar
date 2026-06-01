/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { BoxCage } from '../../mol-geo/primitive/box';
import { Box3D, Sphere3D } from '../../mol-math/geometry';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { shapeFromObj } from '../../mol-model-formats/shape/obj';
import { shapeFromPly } from '../../mol-model-formats/shape/ply';
import { Shape } from '../../mol-model/shape';
import { Task } from '../../mol-task';
import { ColorNames } from '../../mol-util/color/names';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject as SO, PluginStateTransform } from '../objects';

export { BoxShape3D };
type BoxShape3D = typeof BoxShape3D
const BoxShape3D = PluginStateTransform.BuiltIn({
    name: 'box-shape-3d',
    display: 'Box Shape',
    from: SO.Root,
    to: SO.Shape.Provider,
    params: {
        bottomLeft: PD.Vec3(Vec3()),
        topRight: PD.Vec3(Vec3.create(1, 1, 1)),
        radius: PD.Numeric(0.15, { min: 0.01, max: 4, step: 0.01 }),
        color: PD.Color(ColorNames.red)
    }
})({
    canAutoUpdate() {
        return true;
    },
    apply({ params }) {
        return Task.create('Shape Representation', async ctx => {
            return new SO.Shape.Provider({
                label: 'Box',
                data: params,
                params: Mesh.Params,
                getShape: (_, data: typeof params) => {
                    const mesh = getBoxMesh(Box3D.create(params.bottomLeft, params.topRight), params.radius);
                    return Shape.create('Box', data, mesh, () => data.color, () => 1, () => 'Box');
                },
                geometryUtils: Mesh.Utils
            }, { label: 'Box' });
        });
    }
});

export function getBoxMesh(box: Box3D, radius: number, oldMesh?: Mesh) {
    const diag = Vec3.sub(Vec3(), box.max, box.min);
    const translateUnit = Mat4.fromTranslation(Mat4(), Vec3.create(0.5, 0.5, 0.5));
    const scale = Mat4.fromScaling(Mat4(), diag);
    const translate = Mat4.fromTranslation(Mat4(), box.min);
    const transform = Mat4.mul3(Mat4(), translate, scale, translateUnit);

    // TODO: optimize?
    const state = MeshBuilder.createState(256, 128, oldMesh);
    state.currentGroup = 1;
    MeshBuilder.addCage(state, transform, BoxCage(), radius, 2, 20);
    const mesh = MeshBuilder.getMesh(state);

    const center = Vec3.scaleAndAdd(Vec3(), box.min, diag, 0.5);
    const sphereRadius = Vec3.distance(box.min, center);
    mesh.setBoundingSphere(Sphere3D.create(center, sphereRadius));

    return mesh;
}

export { ShapeFromPly };
type ShapeFromPly = typeof ShapeFromPly
const ShapeFromPly = PluginStateTransform.BuiltIn({
    name: 'shape-from-ply',
    display: { name: 'Shape from PLY', description: 'Create Shape from PLY data' },
    from: SO.Format.Ply,
    to: SO.Shape.Provider,
    params(a) {
        return {
            transforms: PD.Optional(PD.Value([Mat4.identity()], { isHidden: true })),
            label: PD.Optional(PD.Text('', { isHidden: true }))
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create shape from PLY', async ctx => {
            const shape = await shapeFromPly(a.data, params).runInContext(ctx);
            const props = { label: params.label || 'Shape' };
            return new SO.Shape.Provider(shape, props);
        });
    }
});

export { ShapeFromObj };
type ShapeFromObj = typeof ShapeFromObj
const ShapeFromObj = PluginStateTransform.BuiltIn({
    name: 'shape-from-obj',
    display: { name: 'Shape from OBJ', description: 'Create Shape from OBJ data' },
    from: SO.Format.Obj,
    to: SO.Shape.Provider,
    params(a) {
        return {
            transforms: PD.Optional(PD.Value([Mat4.identity()], { isHidden: true })),
            label: PD.Optional(PD.Text('', { isHidden: true }))
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create shape from OBJ', async ctx => {
            const shape = await shapeFromObj(a.data, params).runInContext(ctx);
            const props = { label: params.label || 'Shape' };
            return new SO.Shape.Provider(shape, props);
        });
    }
});
