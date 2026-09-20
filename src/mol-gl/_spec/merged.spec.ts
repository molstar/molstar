/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRenderObject, createMergedRenderObject, canMergeRenderObjects } from '../render-object';
import { Scene } from '../scene';
import { getGLContext, tryGetGLContext } from './gl';
import { setDebugMode } from '../../mol-util/debug';
import { ColorNames } from '../../mol-util/color/names';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { ValueCell } from '../../mol-util';

function createSpheresRenderObject(centers: Float32Array, groups: Float32Array, count: number) {
    const spheres = Spheres.create(centers, groups, count);
    const props = PD.getDefaultValues(Spheres.Params);
    const values = Spheres.Utils.createValuesSimple(spheres, props, ColorNames.orange, 1);
    const state = Spheres.Utils.createRenderableState(props);
    return createRenderObject('spheres', values, state, -1);
}

function createMeshRenderObject(vertexCount: number, offset: number) {
    const vertices = new Float32Array(vertexCount * 3);
    const normals = new Float32Array(vertexCount * 3);
    const groups = new Float32Array(vertexCount);
    for (let i = 0; i < vertexCount; ++i) {
        vertices[i * 3] = i + offset;
        normals[i * 3 + 2] = 1;
    }
    const triangleCount = Math.floor(vertexCount / 3);
    const indices = new Uint32Array(triangleCount * 3);
    for (let i = 0, il = triangleCount * 3; i < il; ++i) indices[i] = i;
    const mesh = Mesh.create(vertices, indices, normals, groups, vertexCount, triangleCount);
    const props = PD.getDefaultValues(Mesh.Params);
    const values = Mesh.Utils.createValuesSimple(mesh, props, ColorNames.orange, 1);
    const state = Mesh.Utils.createRenderableState(props);
    return createRenderObject('mesh', values, state, -1);
}

function createMembers() {
    const a = createSpheresRenderObject(new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0]), new Float32Array([0, 1, 2]), 3);
    const b = createSpheresRenderObject(new Float32Array([2, 0, 0, 0, 2, 0]), new Float32Array([0, 1]), 2);
    return [a, b];
}

describe('merged', () => {
    it('values', () => {
        const [a, b] = createMembers();
        expect(canMergeRenderObjects([a, b])).toBe(true);

        const merged = createMergedRenderObject([a, b]);
        const mv = merged.values;

        expect(mv.drawCount.ref.value).toBe(5 * 6);
        expect(mv.uVertexCount.ref.value).toBe(5 * 6);
        expect(mv.instanceCount.ref.value).toBe(2);
        expect(mv.uInstanceCount.ref.value).toBe(2);

        const { segments } = merged.merged;
        expect(segments.count).toBe(2);
        expect(segments.instanceBases).toEqual([0, 1]);
        expect(segments.instanceCounts).toEqual([1, 1]);
        expect(segments.vertexBases).toEqual([0, 18]);

        expect(Array.from(mv.aInstance.ref.value.subarray(0, 2))).toEqual([0, 1]);

        const g0 = a.values.uGroupCount.ref.value;
        const g1 = b.values.uGroupCount.ref.value;
        const aSegment = (mv as any).aSegment.ref.value as Float32Array;
        expect(Array.from(aSegment.subarray(0, 6))).toEqual([
            g0, 0, 0,
            g1, 1 * g0 - 1 * g1, g0,
        ]);

        // second member's positions start at texel 3 in the merged position texture
        const pg = mv.tPositionGroup.ref.value.array;
        expect(Array.from(pg.subarray(3 * 4, 3 * 4 + 4))).toEqual([2, 0, 0, 0]);
        expect(Array.from(pg.subarray(4 * 4, 4 * 4 + 4))).toEqual([0, 2, 0, 1]);

        expect((mv as any).dSegmented.ref.value).toBe(true);
    });

    it('sync', () => {
        const [a, b] = createMembers();
        const merged = createMergedRenderObject([a, b]);

        // member texture update flows into the merged values
        const pgB = b.values.tPositionGroup.ref.value;
        pgB.array.set([3, 0, 0, 0], 0);
        ValueCell.update(b.values.tPositionGroup, pgB);

        merged.merged.sync();
        const pg = merged.values.tPositionGroup.ref.value.array;
        expect(Array.from(pg.subarray(3 * 4, 3 * 4 + 4))).toEqual([3, 0, 0, 0]);
    });

    it('mesh', () => {
        const a = createMeshRenderObject(6, 0);
        const b = createMeshRenderObject(9, 100);
        expect(canMergeRenderObjects([a, b])).toBe(true);

        const merged = createMergedRenderObject([a, b]);
        const mv = merged.values;

        expect(mv.drawCount.ref.value).toBe(6 + 9);
        expect(mv.uVertexCount.ref.value).toBe(6 + 9);
        expect(mv.instanceCount.ref.value).toBe(2);

        const { segments } = merged.merged;
        expect(segments.vertexBases).toEqual([0, 6]);
        expect(segments.elementOffsets).toEqual([0, 6 * 4]);

        // indices of the second member are rewritten to absolute vertex ids
        const elements = mv.elements.ref.value as Uint32Array;
        expect(Array.from(elements.subarray(0, 6))).toEqual([0, 1, 2, 3, 4, 5]);
        expect(Array.from(elements.subarray(6, 15))).toEqual([6, 7, 8, 9, 10, 11, 12, 13, 14]);

        // concatenated vertex positions
        const positions = mv.aPosition.ref.value as Float32Array;
        expect(positions[0]).toBe(0);
        expect(positions[6 * 3]).toBe(100);
    });

    const ctx = tryGetGLContext(32, 32, { fragDepth: true, textureFloat: true });

    (ctx ? it : it.skip)('scene', async () => {
        const ctx = getGLContext(32, 32);
        const scene = Scene.create(ctx);
        const [a, b] = createMembers();
        const merged = createMergedRenderObject([a, b]);
        scene.add(merged);
        const mergedMesh = createMergedRenderObject([createMeshRenderObject(6, 0), createMeshRenderObject(9, 100)]);
        scene.add(mergedMesh);
        setDebugMode(true);
        expect(() => scene.commit()).not.toThrow();
        expect(() => scene.update(undefined, false)).not.toThrow();
        setDebugMode(false);
        ctx.destroy();
    });
});
