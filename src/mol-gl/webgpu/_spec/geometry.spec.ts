/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { TextureMesh } from '../../../mol-geo/geometry/texture-mesh/texture-mesh';
import { Sphere3D } from '../../../mol-math/geometry';
import { WebGPUTextureData } from '../texture-data';
import { PositionLocation } from '../../../mol-geo/util/location-iterator';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Points } from '../../../mol-geo/geometry/points/points';
import { Spheres } from '../../../mol-geo/geometry/spheres/spheres';
import { Cylinders } from '../../../mol-geo/geometry/cylinders/cylinders';
import { CylindersBuilder } from '../../../mol-geo/geometry/cylinders/cylinders-builder';
import { createTransform } from '../../../mol-geo/geometry/transform-data';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Color } from '../../../mol-util/color';
import { packIntToRGBArray } from '../../../mol-util/number-packing';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ValueCell } from '../../../mol-util/value-cell';
import { createRenderObject } from '../../render-object';
import { createWebGPUGeometry, createWebGPUThemes, createWebGPUOverpaintThemes, WebGPUThemeStride } from '../geometry';
import { calcMeshColorSmoothing } from '../../../mol-geo/geometry/mesh/color-smoothing';

function meshObject(instances = 1) {
    const mesh = Mesh.create(new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0]), new Uint32Array([0, 1, 2]), new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), new Float32Array([0, 1, 1]), 3, 1);
    const props = PD.getDefaultValues(Mesh.Params);
    const matrices = new Float32Array(instances * 16);
    for (let i = 0; i < instances; i++) matrices.set(Mat4.identity(), i * 16);
    const values = Mesh.Utils.createValuesSimple(mesh, props, Color(0xff0000), 1, createTransform(matrices, instances, undefined, 0, 0));
    return createRenderObject('mesh', values, Mesh.Utils.createRenderableState(props), -1);
}

describe('native WebGPU geometry', () => {
    it('ignores cached emission data when the overlay is disabled, retaining uniform emission', () => {
        const object = meshObject();
        ValueCell.update(object.values.tEmissive, { array: new Uint8Array([255, 128]), width: 2, height: 1 });
        ValueCell.update(object.values.uEmissive, 0.25); ValueCell.update(object.values.uEmissiveStrength, 0.5);
        const geometry = createWebGPUGeometry(object);
        expect(createWebGPUThemes(object, geometry)[6]).toBe(0.25);
        ValueCell.update(object.values.dEmissive, true);
        expect(createWebGPUThemes(object, geometry)[6]).toBe(0.75);
        ValueCell.update(object.values.dEmissive, false);
        expect(createWebGPUThemes(object, geometry)[6]).toBe(0.25);
    });

    it('keeps group overpaint separate from spatial colors and uses reordered instance IDs', () => {
        const object = meshObject(2);
        ValueCell.update(object.values.dColorType, 'volumeInstance');
        ValueCell.update(object.values.dOverpaint, true); ValueCell.update(object.values.uGroupCount, 2);
        ValueCell.update(object.values.aInstance, new Float32Array([1, 0]));
        ValueCell.update(object.values.tOverpaint, { array: new Uint8Array([255, 0, 0, 128, 0, 255, 0, 255, 0, 0, 255, 64, 255, 255, 0, 0]), width: 4, height: 1 });
        const geometry = createWebGPUGeometry(object);
        const data = createWebGPUOverpaintThemes(object, geometry);
        expect(Array.from(data.subarray(0, 3))).toEqual([0, 0, 1]);
        expect(data[3]).toBeCloseTo(64 / 255);
        expect(Array.from(data.subarray(12, 15))).toEqual([1, 0, 0]);
        expect(data[15]).toBeCloseTo(128 / 255);
        expect(() => createWebGPUThemes(object, geometry)).not.toThrow();
    });

    it.each([1, 3] as const)('packs %i-channel smoothing grids into native RGBA textures', itemSize => {
        const texture = new WebGPUTextureData();
        const result = calcMeshColorSmoothing({
            vertexCount: 1, instanceCount: 1, groupCount: 1,
            transformBuffer: new Float32Array(Mat4.identity()), instanceBuffer: new Float32Array([0]),
            positionBuffer: new Float32Array(3), groupBuffer: new Float32Array(1),
            colorData: { width: 1, height: 1, array: new Uint8Array(itemSize === 1 ? [128] : [10, 20, 30]) },
            colorType: 'groupInstance', boundingSphere: Sphere3D.create(Vec3(), 1),
            invariantBoundingSphere: Sphere3D.create(Vec3(), 1), itemSize
        }, { resolution: 1, stride: 1 }, undefined, texture);
        expect(result.kind).toBe('volume');
        const data = texture.data.array;
        let populated = 0;
        for (let i = 0; i < data.length; i += 4) {
            if (itemSize === 1 ? data[i + 3] > 0 : data[i] > 0) {
                expect(Array.from(data.subarray(i, i + 4))).toEqual(itemSize === 1 ? [0, 0, 0, 128] : [10, 20, 30, 255]);
                populated++;
            }
        }
        expect(populated).toBeGreaterThan(0);
        texture.destroy();
    });

    it('builds spatial material grids using reordered logical instance IDs without WebGL', () => {
        const texture = new WebGPUTextureData(), transforms = new Float32Array(32);
        transforms.set(Mat4.fromTranslation(Mat4(), Vec3.create(-5, 0, 0)));
        transforms.set(Mat4.fromTranslation(Mat4(), Vec3.create(5, 0, 0)), 16);
        const result = calcMeshColorSmoothing({
            vertexCount: 1, instanceCount: 2, groupCount: 1,
            transformBuffer: transforms, instanceBuffer: new Float32Array([1, 0]),
            positionBuffer: new Float32Array(3), groupBuffer: new Float32Array(1),
            colorData: { width: 2, height: 1, array: new Uint8Array([255, 0, 0, 255, 0, 255, 0, 255]) },
            colorType: 'groupInstance', boundingSphere: Sphere3D.create(Vec3(), 10),
            invariantBoundingSphere: Sphere3D.create(Vec3(), 1), itemSize: 4
        }, { resolution: 1, stride: 1 }, undefined, texture);
        expect(result.kind).toBe('volume');
        if (result.kind !== 'volume') throw new Error('Expected a native grid.');
        expect(result.type).toBe('volumeInstance');
        for (const [position, rgba] of [[-5, [0, 255, 0, 255]], [5, [255, 0, 0, 255]]] as const) {
            const x = position - result.gridTransform[0], y = -result.gridTransform[1], z = -result.gridTransform[2];
            const nx = result.gridDim[0], ny = result.gridDim[1], width = texture.data.width;
            const offset = ((Math.floor(z * nx / width) * ny + y) * width + z * nx % width + x) * 4;
            expect(Array.from(texture.data.array.subarray(offset, offset + 4))).toEqual(rgba);
        }
        texture.destroy();
    });

    it('decodes texture meshes and exposes transformed position-theme locations without WebGL', () => {
        const positions = new WebGPUTextureData(), normals = new WebGPUTextureData(), groups = new WebGPUTextureData();
        positions.load({ width: 2, height: 2, array: new Float32Array([0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 99, 99, 99, 1]) });
        normals.load({ width: 2, height: 2, array: new Float32Array([0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0]) });
        const encoded = new Uint8Array(16);
        for (let i = 0; i < 3; i++) packIntToRGBArray(65537 + i, encoded, i * 4);
        groups.load({ width: 2, height: 2, array: encoded });
        const mesh = TextureMesh.create(3, 65540, positions, groups, normals, Sphere3D.create(Vec3(), 2));
        const props = PD.getDefaultValues(TextureMesh.Params);
        const transform = Mat4.fromTranslation(Mat4(), Vec3.create(3, 4, 5));
        const transforms = createTransform(new Float32Array(transform), 1, undefined, 0, 0);
        const object = createRenderObject('texture-mesh', TextureMesh.Utils.createValuesSimple(mesh, props, Color(0xff0000), 1, transforms), TextureMesh.Utils.createRenderableState(props), -1);
        const geometry = createWebGPUGeometry(object);
        expect(geometry.vertexCount).toBe(3);
        expect(Array.from(geometry.indices)).toEqual([0, 1, 2]);
        expect([geometry.vertices[16], geometry.vertices[36], geometry.vertices[56]]).toEqual([65537, 65538, 65539]);
        const iterator = TextureMesh.Utils.createPositionIterator(mesh, transforms);
        expect(iterator.count).toBe(3);
        expect((iterator.move().location as PositionLocation).position).toEqual(Vec3.create(3, 4, 5));
        expect((iterator.move().location as PositionLocation).position).toEqual(Vec3.create(4, 4, 5));
        const floats = Float32Array.from(encoded, v => v / 255);
        groups.load({ width: 2, height: 2, array: floats });
        expect(createWebGPUGeometry(object).vertices).toEqual(geometry.vertices);
        ValueCell.update(object.values.drawCount, 4);
        expect(() => createWebGPUGeometry(object)).toThrow(/complete triangles/);
        positions.destroy(); groups.destroy(); normals.destroy();
    });

    it('preserves mesh attributes and indexes', () => {
        const object = meshObject();
        const geometry = createWebGPUGeometry(object);
        expect(geometry.indices).toEqual(new Uint32Array([0, 1, 2]));
        expect(geometry.vertexCount).toBe(3);
        expect(Array.from(geometry.vertices.subarray(20, 24))).toEqual([1, 0, 0, 1]);
        expect(geometry.vertices[36]).toBe(1);
    });

    it('uses logical instance ids after instance-grid reordering', () => {
        const object = meshObject(2);
        ValueCell.update(object.values.aInstance, new Float32Array([1, 0]));
        ValueCell.update(object.values.uGroupCount, 2);
        ValueCell.update(object.values.dColorType, 'groupInstance');
        ValueCell.update(object.values.tColor, { array: new Uint8Array([255, 0, 0, 0, 255, 0, 0, 0, 255, 255, 255, 0]), width: 4, height: 1 });
        const themes = createWebGPUThemes(object, createWebGPUGeometry(object));
        expect(Array.from(themes.subarray(0, 3))).toEqual([0, 0, 1]);
        expect(themes[7]).toBe(1);
        expect(Array.from(themes.subarray(36, 39))).toEqual([1, 0, 0]);
        expect(themes[43]).toBe(0);
    });

    it('preserves base colors for fragment overpaint and applies per-group transparency', () => {
        const object = meshObject();
        ValueCell.update(object.values.dOverpaint, true);
        ValueCell.update(object.values.tOverpaint, { array: new Uint8Array([0, 0, 255, 255, 0, 255, 0, 255]), width: 2, height: 1 });
        ValueCell.update(object.values.dTransparency, true);
        ValueCell.update(object.values.tTransparency, { array: new Uint8Array([255, 0]), width: 2, height: 1 });
        object.state.alphaFactor = 0.5;
        const geometry = createWebGPUGeometry(object), themes = createWebGPUThemes(object, geometry);
        expect(Array.from(themes.subarray(0, 4))).toEqual([1, 0, 0, 0]);
        expect(Array.from(themes.subarray(12, 16))).toEqual([1, 0, 0, 0.5]);
        const overpaint = createWebGPUOverpaintThemes(object, geometry);
        expect(Array.from(overpaint.subarray(0, 4))).toEqual([0, 0, 1, 1]);
        expect(Array.from(overpaint.subarray(4, 8))).toEqual([0, 1, 0, 1]);
    });

    it('expands point sprites and decodes packed size values', () => {
        const points = Points.create(new Float32Array([0, 0, 0]), new Float32Array([0]), 1);
        const props = PD.getDefaultValues(Points.Params);
        const values = Points.Utils.createValuesSimple(points, props, Color(0xffffff), 1);
        ValueCell.update(values.dSizeType, 'group');
        ValueCell.update(values.uSizeFactor, 2);
        const packed = new Uint8Array(3);
        packIntToRGBArray(250, packed, 0);
        ValueCell.update(values.tSize, { array: packed, width: 1, height: 1 });
        const object = createRenderObject('points', values, Points.Utils.createRenderableState(props), -1);
        const geometry = createWebGPUGeometry(object);
        expect(geometry.vertexCount).toBe(4);
        expect(geometry.indices.length).toBe(6);
        expect(createWebGPUThemes(object, geometry)[4]).toBe(5);
    });

    it('preserves substance groups, reordered instances and material pre-mixing', () => {
        const object = meshObject(2);
        ValueCell.update(object.values.uGroupCount, 2);
        ValueCell.update(object.values.aInstance, new Float32Array([1, 0]));
        ValueCell.update(object.values.dSubstance, true);
        ValueCell.update(object.values.dSubstanceType, 'groupInstance');
        ValueCell.update(object.values.uSubstanceStrength, 0.5);
        ValueCell.update(object.values.tSubstance, { array: new Uint8Array([255, 0, 255, 255, 0, 0, 0, 0, 255, 0, 0, 255, 0, 255, 0, 128]), width: 4, height: 1 });
        const themes = createWebGPUThemes(object, createWebGPUGeometry(object));
        expect(Array.from(themes.subarray(8, 12))).toEqual([0.5, 0, 0, 0.5]);
        expect(themes[23]).toBeCloseTo(128 / 255 * 0.5);
        expect(Array.from(themes.subarray(44, 48))).toEqual([0.5, 0, 0.5, 0.5]);
        expect(themes[59]).toBe(0);
        ValueCell.update(object.values.dSubstanceType, 'instance');
        const instanceThemes = createWebGPUThemes(object, createWebGPUGeometry(object));
        expect(instanceThemes[11]).toBe(0);
        expect(instanceThemes[23]).toBe(0);
        expect(instanceThemes[47]).toBe(0.5);
    });

    it('converts one impostor sphere once, using its actual sphere count', () => {
        const sphere = Spheres.create(new Float32Array([1, 2, 3]), new Float32Array([9]), 1);
        const props = PD.getDefaultValues(Spheres.Params);
        const values = Spheres.Utils.createValuesSimple(sphere, props, Color(0xff0000), 2);
        const geometry = createWebGPUGeometry(createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1));
        expect(geometry.vertexCount).toBe(12);
        expect(geometry.indices.length).toBe(42);
        expect(Array.from(geometry.vertices.subarray(0, 3))).toEqual([1, 2, 3]);
    });

    it('retains prepared sphere subsets, groups and radius scales for each distance level', () => {
        const sphere = Spheres.create(new Float32Array([0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0, 5, 0, 0]), new Float32Array([10, 11, 12, 13, 14, 15]), 6);
        const props = { ...PD.getDefaultValues(Spheres.Params), lodLevels: [
            { minDistance: 0, maxDistance: 25, overlap: 0, stride: 1, scaleBias: 1 },
            { minDistance: 25, maxDistance: 100, overlap: 0, stride: 2, scaleBias: 1 },
        ] };
        const values = Spheres.Utils.createValuesSimple(sphere, props, Color(0xff0000), 1);
        const geometry = createWebGPUGeometry(createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1));
        expect(geometry.sphereLods).toEqual([{ first: 0, count: 252, min: 0, max: 25 }, { first: 252, count: 126, min: 25, max: 100 }]);
        const groups = Array.from({ length: 9 }, (_, i) => geometry.vertices[i * 12 * 20 + 16]);
        expect(groups).toEqual([10, 12, 14, 11, 13, 15, 10, 12, 14]);
        expect(geometry.vertices[7]).toBe(1);
        expect(geometry.vertices[6 * 12 * 20 + 7]).toBe(2);
    });

    it('preserves dual bond endpoint colors, fixed modes, half-bond gradients and caps', () => {
        for (const mode of [0, 0.25, 1, 2, 3]) {
            const builder = CylindersBuilder.create(1, 1);
            builder.add(0, 0, 0, 2, 0, 0, 1, true, true, mode, 0);
            const props = { ...PD.getDefaultValues(Cylinders.Params), colorMode: 'interpolate' as const };
            const values = Cylinders.Utils.createValuesSimple(builder.getCylinders(), props, Color(0xffffff), 1);
            ValueCell.update(values.dColorType, 'group');
            ValueCell.update(values.tColor, { array: new Uint8Array([255, 0, 0, 0, 0, 255]), width: 2, height: 1 });
            const object = createRenderObject('cylinders', values, Cylinders.Utils.createRenderableState(props), -1);
            const geometry = createWebGPUGeometry(object);
            const themes = createWebGPUThemes(object, geometry);
            for (let j = 0; j < geometry.vertexCount; j++) {
                if (geometry.vertices[j * 20 + 15] < 0) continue;
                const weight = mode === 3 ? geometry.vertices[j * 20 + 13] * 0.5 : mode === 2 ? 0 : mode;
                expect(Array.from(themes.subarray(j * WebGPUThemeStride, j * WebGPUThemeStride + 3))).toEqual([1 - weight, 0, weight]);
            }
        }
    });

    it('winds cylinder walls and both end caps toward their outward normals', () => {
        const builder = CylindersBuilder.create(1, 1);
        builder.add(0, 0, 0, 2, 1, 0, 1, true, true, 2, 0);
        const props = PD.getDefaultValues(Cylinders.Params);
        const values = Cylinders.Utils.createValuesSimple(builder.getCylinders(), props, Color(0xffffff), 1);
        const geometry = createWebGPUGeometry(createRenderObject('cylinders', values, Cylinders.Utils.createRenderableState(props), -1));
        const position = (index: number) => Vec3.create(...[0, 1, 2].map(c => geometry.vertices[index * 20 + c] + geometry.vertices[index * 20 + 8 + c]) as [number, number, number]);
        for (let i = 0; i < geometry.indices.length; i += 3) {
            const [a, b, c] = Array.from(geometry.indices.subarray(i, i + 3));
            if (geometry.vertices[a * 20 + 15] < 0) continue; // Camera-plane patches are positioned by the vertex shader.
            const normal = Vec3.fromArray(Vec3(), geometry.vertices, a * 20 + 4);
            const face = Vec3.cross(Vec3(), Vec3.sub(Vec3(), position(b), position(a)), Vec3.sub(Vec3(), position(c), position(a)));
            expect(Vec3.dot(face, normal)).toBeGreaterThan(0);
        }
    });

    it('rejects invalid index counts and unsupported theme formats explicitly', () => {
        const object = meshObject();
        ValueCell.update(object.values.drawCount, 100);
        expect(() => createWebGPUGeometry(object)).toThrow('draw count exceeds');
        ValueCell.update(object.values.drawCount, 3);
        ValueCell.update(object.values.dColorType, 'unsupported');
        ValueCell.update(object.values.uColor, Vec3());
        expect(() => createWebGPUThemes(object, createWebGPUGeometry(object))).toThrow('granularity');
    });
});
