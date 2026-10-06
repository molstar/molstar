/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { Tensor, Mat4 } from '../../../mol-math/linear-algebra';
import { Volume } from '../../../mol-model/volume';
import { CustomProperties } from '../../../mol-model/custom-property';
import { createVolumeTexture3d } from '../../../mol-repr/volume/util';
import { fromHalfFloat } from '../../../mol-util/number-conversion';
import { WebGPUTextureData } from '../texture-data';

function constantVolume(): Volume {
    const space = Tensor.Space([2, 3, 4], [1, 0, 2], Float32Array);
    const data = space.create(); data.fill(7);
    return {
        grid: { cells: Tensor.create(space, data), transform: { kind: 'matrix', matrix: Mat4.identity() }, stats: { min: 7, max: 7, mean: 7, sigma: 0 } },
        instances: [{ transform: Mat4.identity() }], sourceData: { kind: 'test', name: 'constant', data: undefined },
        customProperties: new CustomProperties(), _propertyData: {}, _localPropertyData: {},
    };
}

describe('native volume data', () => {
    it('normalizes constant-density grids without NaN in any upload format', () => {
        const volume = constantVolume();
        for (const format of ['byte', 'float', 'halfFloat'] as const) {
            const image = createVolumeTexture3d(volume, format);
            expect([image.width, image.height, image.depth]).toEqual([2, 3, 4]);
            for (let i = 3; i < image.array.length; i += 4) {
                const density = format === 'halfFloat' ? fromHalfFloat(image.array[i]) : image.array[i];
                expect(density).toBe(0);
            }
        }
        volume.customProperties.dispose();
    });
    it('owns CPU data without allocating a WebGL texture and rejects invalid uploads', () => {
        const texture = new WebGPUTextureData();
        expect(() => texture.load({ width: 2, height: 2, depth: 2, array: new Float32Array(4) })).toThrow(/dimensions/);
        const data = { width: 2, height: 2, depth: 2, array: new Float32Array(32) };
        texture.load(data);
        expect(texture.data).toBe(data);
        expect(texture.getByteCount()).toBe(128);
        expect(() => texture.bind()).toThrow(/WebGL/);
        texture.destroy();
        expect(() => texture.data).toThrow(/unavailable/);
        expect(() => texture.load(data)).toThrow(/destroyed/);
    });
});
