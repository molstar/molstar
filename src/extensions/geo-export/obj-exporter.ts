/**
 * Copyright (c) 2021-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

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
            const [r, g, b] = Color.toRgbNormalized(color).map(v => Math.round(v * 1000) / 1000);
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

    protected async addMeshWithColors(input: AddMeshInput) {
        const { mesh, values, isGeoTexture, mode, webgl, ctx } = input;
        if (mode !== 'triangles') return;

        const obj = this.obj;
        const t = Mat4();
        const n = Mat3();
        const tmpV = Vec3();
        const stride = isGeoTexture ? 4 : 3;

        const colorType = values.dColorType.ref.value;
        const overpaintType = values.dOverpaintType.ref.value;
        const transparencyType = values.dTransparencyType.ref.value;
        const uAlpha = values.uAlpha.ref.value;
        const aTransform = values.aTransform.ref.value;
        const instanceCount = values.uInstanceCount.ref.value;

        let interpolatedColors: Uint8Array | undefined;
        if (webgl && mesh && (colorType === 'volume' || colorType === 'volumeInstance')) {
            interpolatedColors = ObjExporter.getInterpolatedColors(webgl, { vertices: mesh.vertices, vertexCount: mesh.vertexCount, values, stride, colorType });
        }

        let interpolatedOverpaint: Uint8Array | undefined;
        if (webgl && mesh && overpaintType === 'volumeInstance') {
            interpolatedOverpaint = ObjExporter.getInterpolatedOverpaint(webgl, { vertices: mesh.vertices, vertexCount: mesh.vertexCount, values, stride, colorType: overpaintType });
        }

        let interpolatedTransparency: Uint8Array | undefined;
        if (webgl && mesh && transparencyType === 'volumeInstance') {
            const stride = isGeoTexture ? 4 : 3;
            interpolatedTransparency = ObjExporter.getInterpolatedTransparency(webgl, { vertices: mesh.vertices, vertexCount: mesh.vertexCount, values, stride, colorType: transparencyType });
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
                v3transformMat3(tmpV, v3fromArray(tmpV, normals!, i * stride), n);
                StringBuilder.writeSafe(obj, 'vn ');
                StringBuilder.writeFloat(obj, tmpV[0], 100);
                StringBuilder.whitespace1(obj);
                StringBuilder.writeFloat(obj, tmpV[1], 100);
                StringBuilder.whitespace1(obj);
                StringBuilder.writeFloat(obj, tmpV[2], 100);
                StringBuilder.newline(obj);
            }

            const geoData = { values, groups, vertexCount, instanceIndex, isGeoTexture, mode };

            // color
            const quantizedColors = new Uint8Array(drawCount * 3);
            for (let i = 0; i < drawCount; i += 3) {
                const v = isGeoTexture ? i : indices![i];
                const color = ObjExporter.getColor(v, geoData, interpolatedColors, interpolatedOverpaint);
                Color.toArray(color, quantizedColors, i);
            }
            ObjExporter.quantizeColors(quantizedColors, vertexCount);

            // face
            for (let i = 0; i < drawCount; i += 3) {
                const color = Color.fromArray(quantizedColors, i);

                const transparency = ObjExporter.getTransparency(i, geoData, interpolatedTransparency);
                const alpha = Math.round(uAlpha * (1 - transparency) * 10) / 10; // quantized

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

    async getData() {
        return {
            obj: StringBuilder.getString(this.obj),
            mtl: StringBuilder.getString(this.mtl)
        };
    }

    async getBlob(ctx: RuntimeContext) {
        const { obj, mtl } = await this.getData();
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