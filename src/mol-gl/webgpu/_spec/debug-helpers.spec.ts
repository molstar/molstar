/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { BoundingSphereHelper } from '../../../extensions/debug-helpers/bounding-sphere-helper';
import { MeshHelper } from '../../../extensions/debug-helpers/mesh-helper';
import { ClipObjectHelper } from '../../../extensions/debug-helpers/clip-object-helper';
import { DebugHelperParent, DebugRegistry } from '../../../mol-canvas3d/helper/debug-registry';
import { WebGPUHelperScene } from '../../../mol-canvas3d/helper/webgpu-scene';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Sphere3D } from '../../../mol-math/geometry';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Color } from '../../../mol-util/color';
import { Clip } from '../../../mol-util/clip';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ValueCell } from '../../../mol-util/value-cell';
import { createRenderObject, GraphicsRenderObject } from '../../render-object';

function fixture() {
    const mesh = Mesh.create(new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0]), new Uint32Array([0, 1, 2]), new Float32Array([1, 1, 0, 1, 1, 0, 1, 1, 0]), new Float32Array(3), 3, 1);
    const props = PD.getDefaultValues(Mesh.Params);
    const object = createRenderObject('mesh', Mesh.Utils.createValuesSimple(mesh, props, Color(0xffffff), 1), Mesh.Utils.createRenderableState(props), -1);
    const objects: GraphicsRenderObject[] = [object];
    const parent: DebugHelperParent = {
        boundingSphere: Sphere3D.create(Vec3(), 2), boundingSphereVisible: Sphere3D.create(Vec3(), 2),
        has: o => objects.includes(o), forEach: callback => objects.forEach(o => callback({ values: o.values }, o)),
    };
    return { object, objects, parent };
}

describe('CPU-owned native debug helpers', () => {
    it('refreshes normal geometry on transform changes, uses inverse-transpose normals, and removes stale objects', () => {
        const { object, objects, parent } = fixture();
        const scene = new WebGPUHelperScene(), helper = new MeshHelper(undefined, parent, { meshNormals: true }, scene);
        helper.update(); helper.syncVisibility();
        expect(scene.renderObjects).toHaveLength(1);
        const old = scene.renderObjects[0];
        const transform = Mat4.fromScaling(Mat4(), Vec3.create(2, 1, 1)); transform[12] = 5;
        ValueCell.update(object.values.aTransform, new Float32Array(transform));
        helper.update(); helper.syncVisibility();
        expect(scene.renderObjects).toHaveLength(1);
        expect(scene.renderObjects[0].id).not.toEqual(old.id);
        const lines = scene.renderObjects[0];
        if (lines.type !== 'lines') throw new Error('Expected normal lines.');
        const start = lines.values.aStart.ref.value, end = lines.values.aEnd.ref.value;
        expect(start[0]).toBeCloseTo(5);
        expect(end[0] - start[0]).toBeCloseTo(0.1 / Math.sqrt(5), 5);
        expect(end[1] - start[1]).toBeCloseTo(0.2 / Math.sqrt(5), 5);
        objects.length = 0; helper.update();
        expect(scene.renderObjects).toHaveLength(0);
    });

    it('rebuilds every sphere category after clear and refreshes sphere geometry bounds', () => {
        const { object, parent } = fixture();
        const scene = new WebGPUHelperScene(), helper = new BoundingSphereHelper(undefined, parent, { sceneBoundingSpheres: true, visibleSceneBoundingSpheres: true, objectBoundingSpheres: true, instanceBoundingSpheres: true }, scene);
        helper.update(); helper.syncVisibility();
        expect(scene.renderObjects).toHaveLength(4);
        const sphere = scene.renderObjects[0];
        const oldRadius = sphere.values.invariantBoundingSphere.ref.value.radius;
        parent.boundingSphere.radius = 4; helper.update();
        expect(sphere.values.invariantBoundingSphere.ref.value.radius).toBeCloseTo(oldRadius * 2);
        helper.clear(); expect(scene.renderObjects).toHaveLength(0);
        helper.update(); helper.syncVisibility();
        expect(scene.renderObjects).toHaveLength(4);
        expect(scene.renderObjects.every(o => o.state.visible)).toBe(true);
        object.state.visible = false; helper.syncVisibility();
        expect(scene.renderObjects.filter(o => o.state.visible)).toHaveLength(2);
    });

    it('clears replaced and unregistered scenes without retaining geometry', () => {
        const { parent } = fixture();
        const registry = new DebugRegistry<WebGPUHelperScene, undefined, DebugHelperParent>(undefined, parent);
        const first = new MeshHelper(undefined, parent, { meshNormals: true }, new WebGPUHelperScene());
        const second = new MeshHelper(undefined, parent, { meshNormals: true }, new WebGPUHelperScene());
        registry.register('mesh', first); registry.update();
        expect(first.scene.renderObjects).toHaveLength(1);
        registry.register('mesh', second); registry.update();
        expect(first.scene.renderObjects).toHaveLength(0);
        expect(second.scene.renderObjects).toHaveLength(1);
        registry.unregister('mesh');
        expect(second.scene.renderObjects).toHaveLength(0);
        expect(registry.scenes).toHaveLength(0);
        expect(registry.isEnabled).toBe(false);
    });

    it('resizes unbounded clip helpers when the parent scene grows', () => {
        for (const type of ['plane', 'infiniteCone'] as const) {
            const { object, parent } = fixture();
            const props = { ...PD.getDefaultValues(Mesh.Params), clip: { variant: 'pixel' as const, objects: [{ ...Clip.Params.objects.ctor(), type, scale: Vec3.create(1, 1, 1) }] } };
            Mesh.Utils.updateValues(object.values, props);
            parent.boundingSphereVisible.radius = 20;
            const scene = new WebGPUHelperScene(), helper = new ClipObjectHelper(undefined, parent, { clipObjects: true }, scene);
            helper.update();
            expect(scene.renderObjects).toHaveLength(2);
            const original = scene.renderObjects[0];
            parent.boundingSphereVisible.radius = 40; helper.update();
            expect(scene.renderObjects).toHaveLength(2);
            const updated = scene.renderObjects[0];
            expect(updated.id).not.toBe(original.id);
            expect(updated.values.invariantBoundingSphere.ref.value.radius).toBeCloseTo(original.values.invariantBoundingSphere.ref.value.radius * 2);
        }
    });
});
