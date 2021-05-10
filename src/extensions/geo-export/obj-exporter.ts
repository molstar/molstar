/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import { BaseValues } from '../../mol-gl/renderable/schema';
import { asciiWrite } from '../../mol-io/common/ascii';
import { Vec3, Mat3, Mat4 } from '../../mol-math/linear-algebra';
import { RuntimeContext } from '../../mol-task';
import { StringBuilder } from '../../mol-util';
import { Color } from '../../mol-util/color/color';
import { zip } from '../../mol-util/zip/zip';
import { MeshExporter } from './mesh-exporter';

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

    protected async addMeshWithColors(vertices: Float32Array, normals: Float32Array, indices: Uint32Array | undefined, groups: Float32Array | Uint8Array, vertexCount: number, drawCount: number, values: BaseValues, instanceIndex: number, isGeoTexture: boolean, ctx: RuntimeContext) {
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

        Mat4.fromArray(t, aTransform, instanceIndex * 16);
        mat3directionTransform(n, t);

        const currentProgress = (vertexCount * 2 + drawCount) * instanceIndex;
        await ctx.update({ isIndeterminate: false, current: currentProgress, max: (vertexCount * 2 + drawCount) * values.uInstanceCount.ref.value });

        // position
        for (let i = 0; i < vertexCount; ++i) {
            if (i % 1000 === 0 && ctx.shouldUpdate) await ctx.update({ current: currentProgress + i });
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
            if (i % 1000 === 0 && ctx.shouldUpdate) await ctx.update({ current: currentProgress + vertexCount + i });
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
            if (i % 3000 === 0 && ctx.shouldUpdate) await ctx.update({ current: currentProgress + vertexCount * 2 + i });
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
                    color = Color.fromArray(tColor, (instanceIndex * drawCount + indices![i]) * 3);
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

    constructor(private filename: string) {
        super();
        StringBuilder.writeSafe(this.obj, `mtllib ${filename}.mtl\n`);
    }
}