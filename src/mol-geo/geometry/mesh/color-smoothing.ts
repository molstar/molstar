/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createTextureImage, TextureImage } from '../../../mol-gl/renderable/util';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { Texture } from '../../../mol-gl/webgl/texture';
import { Box3D, Sphere3D } from '../../../mol-math/geometry';
import { Vec2, Vec3, Vec4 } from '../../../mol-math/linear-algebra';
import { getVolumeTexture2dLayout } from '../../../mol-repr/volume/util';
import { Color } from '../../../mol-util/color';

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

export function calcMeshColorSmoothing(input: ColorSmoothingInput, resolution: number, stride: number, webgl?: WebGLContext, texture?: Texture) {
    const { colorType, vertexCount, groupCount, positionBuffer, transformBuffer, groupBuffer } = input;

    const isInstanceType = colorType.endsWith('Instance');
    const box = Box3D.fromSphere3D(Box3D(), isInstanceType ? input.boundingSphere : input.invariantBoundingSphere);

    const scaleFactor = 1 / resolution;
    const scaledBox = Box3D.scale(Box3D(), box, scaleFactor);
    const gridDim = Box3D.size(Vec3(), scaledBox);
    Vec3.ceil(gridDim, gridDim);
    Vec3.add(gridDim, gridDim, Vec3.create(2, 2, 2));
    const { min } = box;

    const [xn, yn] = gridDim;
    const { width, height } = getVolumeTexture2dLayout(gridDim);
    // console.log({ width, height, dim });

    const itemSize = 3;
    const data = new Float32Array(width * height * itemSize);
    const count = new Float32Array(width * height);

    const grid = new Uint8Array(width * height * itemSize);
    const textureImage: TextureImage<Uint8Array> = { array: grid, width, height, filter: 'linear' };

    const instanceCount = isInstanceType ? input.instanceCount : 1;
    const colors = input.colorData.array;

    function getIndex(x: number, y: number, z: number) {
        const column = Math.floor(((z * xn) % width) / xn);
        const row = Math.floor((z * xn) / width);
        const px = column * xn + x;
        return itemSize * ((row * yn * width) + (y * width) + px);
    }

    const p = 2;
    const [dimX, dimY, dimZ] = gridDim;
    const v = Vec3();

    for (let i = 0; i < instanceCount; ++i) {
        for (let j = 0; j < vertexCount; j += stride) {
            Vec3.fromArray(v, positionBuffer, j * 3);
            if (isInstanceType) Vec3.transformMat4Offset(v, v, transformBuffer, 0, 0, i * 16);
            Vec3.sub(v, v, min);
            Vec3.scale(v, v, scaleFactor);
            const [vx, vy, vz] = v;

            // vertex mapped to grid
            const x = Math.floor(vx);
            const y = Math.floor(vy);
            const z = Math.floor(vz);

            // group colors
            const ci = i * groupCount + groupBuffer[j];
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

                        let s = p - d;
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
        grid[i3] = Math.round(data[i3] / c);
        grid[i3 + 1] = Math.round(data[i3 + 1] / c);
        grid[i3 + 2] = Math.round(data[i3 + 2] / c);
    }

    const gridTexDim = Vec2.create(width, height);
    const gridTransform = Vec4.create(min[0], min[1], min[2], scaleFactor);
    const type = isInstanceType ? 'volumeInstance' as const : 'volume' as const;

    if (webgl) {
        if (!texture) texture = webgl.resources.texture('image-uint8', 'rgb', 'ubyte', 'linear');
        texture.load(textureImage);

        return { kind: 'volume' as const, texture, gridTexDim, gridDim, gridTransform, type };
    } else {
        const interpolated = getTrilinearlyInterpolated({ vertexCount, instanceCount, transformBuffer, positionBuffer, colorType: type, grid, gridDim, gridTexDim, gridTransform, vertexStride: 3, colorStride: 3 });

        return {
            kind: 'vertex' as const,
            texture: interpolated,
            texDim: Vec2.create(interpolated.width, interpolated.height),
            type: isInstanceType ? 'vertexInstance' : 'vertex'
        };
    }
}

//

interface ColorInterpolationInput {
    vertexCount: number
    instanceCount: number
    transformBuffer: Float32Array
    positionBuffer: Float32Array
    colorType: 'volumeInstance' | 'volume'
    grid: Uint8Array // 2d layout
    gridTexDim: Vec2
    gridDim: Vec3
    gridTransform: Vec4
    vertexStride: number
    colorStride: number
}

export function getTrilinearlyInterpolated(input: ColorInterpolationInput): TextureImage<Uint8Array> {
    const { vertexCount, positionBuffer, transformBuffer, grid, gridDim, gridTexDim, gridTransform, vertexStride, colorStride } = input;

    const isInstanceType = input.colorType.endsWith('Instance');
    const instanceCount = isInstanceType ? input.instanceCount : 1;
    const image = createTextureImage(Math.max(1, instanceCount * vertexCount), 3, Uint8Array);
    const { array } = image;

    const [xn, yn] = gridDim;
    const width = gridTexDim[0];
    const min = Vec3.fromArray(Vec3(), gridTransform, 0);
    const scaleFactor = gridTransform[3];

    function getIndex(x: number, y: number, z: number) {
        const column = Math.floor(((z * xn) % width) / xn);
        const row = Math.floor((z * xn) / width);
        const px = column * xn + x;
        return colorStride * ((row * yn * width) + (y * width) + px);
    }

    const v = Vec3();
    const v0 = Vec3();
    const v1 = Vec3();
    const vd = Vec3();

    for (let i = 0; i < instanceCount; ++i) {
        for (let j = 0; j < vertexCount; ++j) {
            Vec3.fromArray(v, positionBuffer, j * vertexStride);
            if (isInstanceType) Vec3.transformMat4Offset(v, v, transformBuffer, 0, 0, i * 16);
            Vec3.sub(v, v, min);
            Vec3.scale(v, v, scaleFactor);

            Vec3.floor(v0, v);
            Vec3.ceil(v1, v);

            Vec3.sub(vd, v, v0);
            Vec3.sub(v, v1, v0);
            Vec3.div(vd, vd, v);

            const [x0, y0, z0] = v0;
            const [x1, y1, z1] = v1;
            const [xd, yd, zd] = vd;

            const s000 = Color.fromArray(grid, getIndex(x0, y0, z0));
            const s100 = Color.fromArray(grid, getIndex(x1, y0, z0));
            const s001 = Color.fromArray(grid, getIndex(x0, y0, z1));
            const s101 = Color.fromArray(grid, getIndex(x1, y0, z1));
            const s010 = Color.fromArray(grid, getIndex(x0, y1, z0));
            const s110 = Color.fromArray(grid, getIndex(x1, y1, z0));
            const s011 = Color.fromArray(grid, getIndex(x0, y1, z1));
            const s111 = Color.fromArray(grid, getIndex(x1, y1, z1));

            const s00 = Color.interpolate(s000, s100, xd);
            const s01 = Color.interpolate(s001, s101, xd);
            const s10 = Color.interpolate(s010, s110, xd);
            const s11 = Color.interpolate(s011, s111, xd);

            const s0 = Color.interpolate(s00, s10, yd);
            const s1 = Color.interpolate(s01, s11, yd);

            Color.toArray(Color.interpolate(s0, s1, zd), array, (i * vertexCount + j) * 3);
        }
    }

    return image;
}
