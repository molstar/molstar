/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import { sort, arraySwap } from '../../mol-data/util';
import { asciiWrite } from '../../mol-io/common/ascii';
import { Box3D } from '../../mol-math/geometry';
import { Vec3, Mat3, Mat4 } from '../../mol-math/linear-algebra';
import { RuntimeContext } from '../../mol-task';
import { StringBuilder } from '../../mol-util';
import { Color } from '../../mol-util/color/color';
import { zip } from '../../mol-util/zip/zip';
import { MeshExporter, AddMeshInput } from './mesh-exporter';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3transformMat4 = Vec3.transformMat4;
const v3transformMat3 = Vec3.transformMat3;
const mat3directionTransform = Mat3.directionTransform;

// http://paulbourke.net/dataformats/obj/
// http://paulbourke.net/dataformats/mtl/

export type ObjData = {
    obj: string
    mtl: string
}

export class ObjExporter extends MeshExporter<ObjData> {
    readonly fileExtension = 'zip';
    private obj = StringBuilder.create();
    private mtl = StringBuilder.create();
    private vertexOffset = 0;
    private currentColor: Color | undefined;
    private currentAlpha: number | undefined;
    private materialSet = new Set<string>();
    private centerTransform: Mat4;

    private updateMaterial(color: Color, alpha: number) {
        if (this.currentColor === color && this.currentAlpha === alpha) return;

        this.currentColor = color;
        this.currentAlpha = alpha;
        const material = Color.toHexString(color) + alpha;
        StringBuilder.writeSafe(this.obj, `usemtl ${material}`);
        StringBuilder.newline(this.obj);
        if (!this.materialSet.has(material)) {
            this.materialSet.add(material);
            const [r, g, b] = Color.toRgbNormalized(color);
            const mtl = this.mtl;
            StringBuilder.writeSafe(mtl, `newmtl ${material}\n`);
            StringBuilder.writeSafe(mtl, 'illum 2\n'); // illumination model
            StringBuilder.writeSafe(mtl, 'Ns 163\n'); // specular exponent
            StringBuilder.writeSafe(mtl, 'Ni 0.001\n'); // optical density a.k.a. index of refraction
            StringBuilder.writeSafe(mtl, 'Ka 0 0 0\n'); // ambient reflectivity
            StringBuilder.writeSafe(mtl, 'Kd '); // diffuse reflectivity
            StringBuilder.writeFloat(mtl, r, 1000);
            StringBuilder.whitespace1(mtl);
            StringBuilder.writeFloat(mtl, g, 1000);
            StringBuilder.whitespace1(mtl);
            StringBuilder.writeFloat(mtl, b, 1000);
            StringBuilder.newline(mtl);
            StringBuilder.writeSafe(mtl, 'Ks 0.25 0.25 0.25\n'); // specular reflectivity
            StringBuilder.writeSafe(mtl, 'd '); // dissolve
            StringBuilder.writeFloat(mtl, alpha, 1000);
            StringBuilder.newline(mtl);
        }
    }

    private static quantizeColors(colorArray: Uint8Array, vertexCount: number) {
        if (vertexCount <= 1024) return;
        const rgb = Vec3();
        const min = Vec3();
        const max = Vec3();
        const sum = Vec3();
        const colorMap = new Map<Color, Color>();
        const colorComparers = [
            (colors: Color[], i: number, j: number) => (Color.toVec3(rgb, colors[i])[0] - Color.toVec3(rgb, colors[j])[0]),
            (colors: Color[], i: number, j: number) => (Color.toVec3(rgb, colors[i])[1] - Color.toVec3(rgb, colors[j])[1]),
            (colors: Color[], i: number, j: number) => (Color.toVec3(rgb, colors[i])[2] - Color.toVec3(rgb, colors[j])[2]),
        ];

        const medianCut = (colors: Color[], l: number, r: number, depth: number) => {
            if (l > r) return;
            if (l === r || depth >= 10) {
                // Find the average color.
                Vec3.set(sum, 0, 0, 0);
                for (let i = l; i <= r; ++i) {
                    Color.toVec3(rgb, colors[i]);
                    Vec3.add(sum, sum, rgb);
                }
                Vec3.round(rgb, Vec3.scale(rgb, sum, 1 / (r - l + 1)));
                const averageColor = Color.fromArray(rgb, 0);
                for (let i = l; i <= r; ++i) colorMap.set(colors[i], averageColor);
                return;
            }

            // Find the color channel with the greatest range.
            Vec3.set(min, 255, 255, 255);
            Vec3.set(max, 0, 0, 0);
            for (let i = l; i <= r; ++i) {
                Color.toVec3(rgb, colors[i]);
                for (let j = 0; j < 3; ++j) {
                    Vec3.min(min, min, rgb);
                    Vec3.max(max, max, rgb);
                }
            }
            let k = 0;
            if (max[1] - min[1] > max[k] - min[k]) k = 1;
            if (max[2] - min[2] > max[k] - min[k]) k = 2;

            sort(colors, l, r + 1, colorComparers[k], arraySwap);

            const m = (l + r) >> 1;
            medianCut(colors, l, m, depth + 1);
            medianCut(colors, m + 1, r, depth + 1);
        };

        // Create an array of unique colors and use the median cut algorithm.
        const colorSet = new Set<Color>();
        for (let i = 0; i < vertexCount; ++i) {
            colorSet.add(Color.fromArray(colorArray, i * 3));
        }
        const colors = Array.from(colorSet);
        medianCut(colors, 0, colors.length - 1, 0);

        // Map actual colors to quantized colors.
        for (let i = 0; i < vertexCount; ++i) {
            const color = colorMap.get(Color.fromArray(colorArray, i * 3));
            Color.toArray(color!, colorArray, i * 3);
        }
    }

    protected async addMeshWithColors(input: AddMeshInput) {
        const { mesh, values, isGeoTexture, webgl, ctx } = input;

        const obj = this.obj;
        const t = Mat4();
        const n = Mat3();
        const tmpV = Vec3();
        const stride = isGeoTexture ? 4 : 3;

        const groupCount = values.uGroupCount.ref.value;
        const colorType = values.dColorType.ref.value;
        const tColor = values.tColor.ref.value.array;
        const uAlpha = values.uAlpha.ref.value;
        const dTransparency = values.dTransparency.ref.value;
        const tTransparency = values.tTransparency.ref.value;
        const aTransform = values.aTransform.ref.value;
        const instanceCount = values.uInstanceCount.ref.value;

        let interpolatedColors: Uint8Array;
        if (colorType === 'volume' || colorType === 'volumeInstance') {
            interpolatedColors = ObjExporter.getInterpolatedColors(mesh!.vertices, mesh!.vertexCount, values, stride, colorType, webgl!);
            ObjExporter.quantizeColors(interpolatedColors, mesh!.vertexCount);
        }

        await ctx.update({ isIndeterminate: false, current: 0, max: instanceCount });

        for (let instanceIndex = 0; instanceIndex < instanceCount; ++instanceIndex) {
            if (ctx.shouldUpdate) await ctx.update({ current: instanceIndex + 1 });

            const { vertices, normals, indices, groups, vertexCount, drawCount } = ObjExporter.getInstance(input, instanceIndex);

            Mat4.fromArray(t, aTransform, instanceIndex * 16);
            Mat4.mul(t, this.centerTransform, t);
            mat3directionTransform(n, t);

            // position
            for (let i = 0; i < vertexCount; ++i) {
                v3transformMat4(tmpV, v3fromArray(tmpV, vertices, i * stride), t);
                StringBuilder.writeSafe(obj, 'v ');
                StringBuilder.writeFloat(obj, tmpV[0], 1000);
                StringBuilder.whitespace1(obj);
                StringBuilder.writeFloat(obj, tmpV[1], 1000);
                StringBuilder.whitespace1(obj);
                StringBuilder.writeFloat(obj, tmpV[2], 1000);
                StringBuilder.newline(obj);
            }

            // normal
            for (let i = 0; i < vertexCount; ++i) {
                v3transformMat3(tmpV, v3fromArray(tmpV, normals, i * stride), n);
                StringBuilder.writeSafe(obj, 'vn ');
                StringBuilder.writeFloat(obj, tmpV[0], 100);
                StringBuilder.whitespace1(obj);
                StringBuilder.writeFloat(obj, tmpV[1], 100);
                StringBuilder.whitespace1(obj);
                StringBuilder.writeFloat(obj, tmpV[2], 100);
                StringBuilder.newline(obj);
            }

            // face
            for (let i = 0; i < drawCount; i += 3) {
                let color: Color;
                switch (colorType) {
                    case 'uniform':
                        color = Color.fromNormalizedArray(values.uColor.ref.value, 0);
                        break;
                    case 'instance':
                        color = Color.fromArray(tColor, instanceIndex * 3);
                        break;
                    case 'group': {
                        const group = isGeoTexture ? ObjExporter.getGroup(groups, i) : groups[indices![i]];
                        color = Color.fromArray(tColor, group * 3);
                        break;
                    }
                    case 'groupInstance': {
                        const group = isGeoTexture ? ObjExporter.getGroup(groups, i) : groups[indices![i]];
                        color = Color.fromArray(tColor, (instanceIndex * groupCount + group) * 3);
                        break;
                    }
                    case 'vertex':
                        color = Color.fromArray(tColor, indices![i] * 3);
                        break;
                    case 'vertexInstance':
                        color = Color.fromArray(tColor, (instanceIndex * vertexCount + indices![i]) * 3);
                        break;
                    case 'volume':
                        color = Color.fromArray(interpolatedColors!, (isGeoTexture ? i : indices![i]) * 3);
                        break;
                    case 'volumeInstance':
                        color = Color.fromArray(interpolatedColors!, (instanceIndex * vertexCount + (isGeoTexture ? i : indices![i])) * 3);
                        break;
                    default: throw new Error('Unsupported color type.');
                }

                let alpha = uAlpha;
                if (dTransparency) {
                    const group = isGeoTexture ? ObjExporter.getGroup(groups, i) : groups[indices![i]];
                    const transparency = tTransparency.array[instanceIndex * groupCount + group] / 255;
                    alpha *= 1 - transparency;
                }

                this.updateMaterial(color, alpha);

                const v1 = this.vertexOffset + (isGeoTexture ? i : indices![i]) + 1;
                const v2 = this.vertexOffset + (isGeoTexture ? i + 1 : indices![i + 1]) + 1;
                const v3 = this.vertexOffset + (isGeoTexture ? i + 2 : indices![i + 2]) + 1;
                StringBuilder.writeSafe(obj, 'f ');
                StringBuilder.writeInteger(obj, v1);
                StringBuilder.writeSafe(obj, '//');
                StringBuilder.writeIntegerAndSpace(obj, v1);
                StringBuilder.writeInteger(obj, v2);
                StringBuilder.writeSafe(obj, '//');
                StringBuilder.writeIntegerAndSpace(obj, v2);
                StringBuilder.writeInteger(obj, v3);
                StringBuilder.writeSafe(obj, '//');
                StringBuilder.writeInteger(obj, v3);
                StringBuilder.newline(obj);
            }

            this.vertexOffset += vertexCount;
        }
    }

    getData() {
        return {
            obj: StringBuilder.getString(this.obj),
            mtl: StringBuilder.getString(this.mtl)
        };
    }

    async getBlob(ctx: RuntimeContext) {
        const { obj, mtl } = this.getData();
        const objData = new Uint8Array(obj.length);
        asciiWrite(objData, obj);
        const mtlData = new Uint8Array(mtl.length);
        asciiWrite(mtlData, mtl);
        const zipDataObj = {
            [this.filename + '.obj']: objData,
            [this.filename + '.mtl']: mtlData
        };
        return new Blob([await zip(ctx, zipDataObj)], { type: 'application/zip' });
    }

    constructor(private filename: string, boundingBox: Box3D) {
        super();
        StringBuilder.writeSafe(this.obj, `mtllib ${filename}.mtl\n`);
        const tmpV = Vec3();
        Vec3.add(tmpV, boundingBox.min, boundingBox.max);
        Vec3.scale(tmpV, tmpV, -0.5);
        this.centerTransform = Mat4.fromTranslation(Mat4(), tmpV);
    }
}