/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { TextureImage } from '../../../mol-gl/renderable/util';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { Box3D, Sphere3D } from '../../../mol-math/geometry';
import { Vec2, Vec3, Vec4 } from '../../../mol-math/linear-algebra';
import { getVolumeTexture2dLayout } from '../../../mol-repr/volume/util';

interface ColorSmoothingInput {
    vertexCount: number
    instanceCount: number
    groupCount: number
    transformBuffer: Float32Array
    instanceBuffer: Float32Array
    positionBuffer: Float32Array
    groupBuffer: Float32Array
    colorData: TextureImage<Uint8Array>
    colorType: 'group' | 'groupInstance'
    boundingSphere: Sphere3D
    invariantBoundingSphere: Sphere3D
}

export function calcMeshColorSmoothing(webgl: WebGLContext, input: ColorSmoothingInput, resolution: number, stride: number) {
    const isInstanceType = input.colorType.endsWith('Instance');
    const box = Box3D.fromSphere3D(Box3D(), isInstanceType ? input.boundingSphere : input.invariantBoundingSphere);

    const scaleFactor = 1 / resolution;
    const scaledBox = Box3D.scale(Box3D(), box, scaleFactor);
    const dim = Box3D.size(Vec3(), scaledBox);
    Vec3.ceil(dim, dim);
    Vec3.add(dim, dim, Vec3.create(2, 2, 2));
    const { min } = box;

    const [ xn, yn ] = dim;
    const { width, height } = getVolumeTexture2dLayout(dim);
    // console.log({ width, height, dim });

    const itemSize = 3;
    const data = new Float32Array(width * height * itemSize);
    const count = new Float32Array(width * height);

    const array = new Uint8Array(width * height * itemSize);
    const textureImage: TextureImage<Uint8Array> = { array, width, height, filter: 'linear' };

    const vertexCount = input.vertexCount;
    const groupCount = input.groupCount;
    const instanceCount = isInstanceType ? input.instanceCount : 1;
    const positions = input.positionBuffer;
    const colors = input.colorData.array;
    const groups = input.groupBuffer;
    const transforms = input.transformBuffer;

    function getIndex(x: number, y: number, z: number) {
        const column = Math.floor(((z * xn) % width) / xn);
        const row = Math.floor((z * xn) / width);
        const px = column * xn + x;
        return itemSize * ((row * yn * width) + (y * width) + px);
    }

    const p = 2;
    const [dimX, dimY, dimZ] = dim;
    const v = Vec3();

    for (let i = 0; i < instanceCount; ++i) {
        for (let j = 0; j < vertexCount; j += stride) {
            Vec3.fromArray(v, positions, j * 3);
            if (isInstanceType) Vec3.transformMat4Offset(v, v, transforms, 0, 0, i * 16);
            Vec3.sub(v, v, min);
            Vec3.scale(v, v, scaleFactor);
            const [vx, vy, vz] = v;

            // vertex mapped to grid
            const x = Math.floor(vx);
            const y = Math.floor(vy);
            const z = Math.floor(vz);

            // group colors
            const ci = i * groupCount + groups[j];
            const r = colors[ci * 3];
            const g = colors[ci * 3 + 1];
            const b = colors[ci * 3 + 2];

            // Extents of grid to consider for this atom
            const begX = Math.max(0, x - p);
            const begY = Math.max(0, y - p);
            const begZ = Math.max(0, z - p);

            // Add two to these points:
            // - x, y, z are floor'd values so this ensures coverage
            // - these are loop limits (exclusive)
            const endX = Math.min(dimX, x + p + 2);
            const endY = Math.min(dimY, y + p + 2);
            const endZ = Math.min(dimZ, z + p + 2);

            for (let xi = begX; xi < endX; ++xi) {
                const dx = xi - vx;
                for (let yi = begY; yi < endY; ++yi) {
                    const dy = yi - vy;
                    for (let zi = begZ; zi < endZ; ++zi) {
                        const dz = zi - vz;
                        const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
                        if (d > p) continue;

                        // const s = Math.pow(p, 10) - Math.pow(d, 10);
                        let s = p - d;
                        // if (d < 0.5) s *= 50;
                        const index = getIndex(xi, yi, zi);
                        data[index] += r * s;
                        data[index + 1] += g * s;
                        data[index + 2] += b * s;
                        count[index / 3] += s;
                    }
                }
            }
        }
    }

    for (let i = 0, il = count.length; i < il; ++i) {
        const i3 = i * 3;
        const c = count[i];
        array[i3] = Math.round(data[i3] / c);
        array[i3 + 1] = Math.round(data[i3 + 1] / c);
        array[i3 + 2] = Math.round(data[i3 + 2] / c);
    }

    const texture = webgl.resources.texture('image-uint8', 'rgb', 'ubyte', 'linear');
    texture.load(textureImage);

    // ValueCell.updateIfChanged(values.dColorType, isInstanceType ? 'volumeInstance' : 'volume');
    // ValueCell.update(values.tColorGrid, texture);
    // ValueCell.update(values.uColorTexDim, Vec2.create(width, height));
    // ValueCell.update(values.uColorGridDim, dim);
    // ValueCell.update(values.uColorGridTransform, Vec4.create(min[0], min[1], min[2], scaleFactor));
    // ValueCell.updateIfChanged(values.dColorGridType, '2d');

    // return values;

    const transform = Vec4.create(min[0], min[1], min[2], scaleFactor);
    const type = isInstanceType ? 'volumeInstance' : 'volume';

    return { texture, gridDim: dim, gridTexDim: Vec2.create(width, height), transform, type };
}