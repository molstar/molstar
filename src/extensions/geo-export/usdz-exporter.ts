/**
 * Copyright (c) 2021-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { asciiWrite } from '../../mol-io/common/ascii';
import { Box3D } from '../../mol-math/geometry';
import { Vec3, Mat3, Mat4 } from '../../mol-math/linear-algebra';
import { PLUGIN_VERSION } from '../../mol-plugin/version';
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

// https://graphics.pixar.com/usd/docs/index.html

export type UsdzData = {
    usdz: ArrayBuffer
}

export class UsdzExporter extends MeshExporter<UsdzData> {
    readonly fileExtension = 'usdz';
    private meshes: string[] = [];
    private materials: string[] = [];
    private materialMap = new Map<string, number>();
    private centerTransform: Mat4;

    private addMaterial(color: Color, alpha: number, metalness: number, roughness: number): number {
        const hash = `${color}|${alpha}|${metalness}|${roughness}`;
        if (this.materialMap.has(hash)) return this.materialMap.get(hash)!;

        const materialKey = this.materialMap.size;
        this.materialMap.set(hash, materialKey);

        const [r, g, b] = Color.toRgbNormalized(color).map(v => Math.round(v * 1000) / 1000);
        this.materials.push(`
def Material "material${materialKey}"
{
    token outputs:surface.connect = </material${materialKey}/shader.outputs:surface>
    def Shader "shader"
    {
        uniform token info:id = "UsdPreviewSurface"
        color3f inputs:diffuseColor = (${r},${g},${b})
        float inputs:opacity = ${alpha}
        float inputs:metallic = ${metalness}
        float inputs:roughness = ${roughness}
        token outputs:surface
    }
}
`);
        return materialKey;
    }

    protected async addMeshWithColors(input: AddMeshInput) {
        const { mesh, values, isGeoTexture, mode, webgl, ctx } = input;
        if (mode !== 'triangles') return;

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
        const metalness = values.uMetalness.ref.value;
        const roughness = values.uRoughness.ref.value;

        let interpolatedColors: Uint8Array | undefined;
        if (webgl && mesh && (colorType === 'volume' || colorType === 'volumeInstance')) {
            interpolatedColors = UsdzExporter.getInterpolatedColors(webgl, { vertices: mesh.vertices, vertexCount: mesh.vertexCount, values, stride, colorType });
        }

        let interpolatedOverpaint: Uint8Array | undefined;
        if (webgl && mesh && overpaintType === 'volumeInstance') {
            const stride = isGeoTexture ? 4 : 3;
            interpolatedOverpaint = UsdzExporter.getInterpolatedOverpaint(webgl, { vertices: mesh.vertices, vertexCount: mesh.vertexCount, values, stride, colorType: overpaintType });
        }

        let interpolatedTransparency: Uint8Array | undefined;
        if (webgl && mesh && transparencyType === 'volumeInstance') {
            const stride = isGeoTexture ? 4 : 3;
            interpolatedTransparency = UsdzExporter.getInterpolatedTransparency(webgl, { vertices: mesh.vertices, vertexCount: mesh.vertexCount, values, stride, colorType: transparencyType });
        }

        await ctx.update({ isIndeterminate: false, current: 0, max: instanceCount });

        for (let instanceIndex = 0; instanceIndex < instanceCount; ++instanceIndex) {
            if (ctx.shouldUpdate) await ctx.update({ current: instanceIndex + 1 });

            const { vertices, normals, indices, groups, vertexCount, drawCount } = UsdzExporter.getInstance(input, instanceIndex);

            Mat4.fromArray(t, aTransform, instanceIndex * 16);
            Mat4.mul(t, this.centerTransform, t);
            mat3directionTransform(n, t);

            const vertexBuilder = StringBuilder.create();
            const normalBuilder = StringBuilder.create();
            const indexBuilder = StringBuilder.create();

            // position
            for (let i = 0; i < vertexCount; ++i) {
                v3transformMat4(tmpV, v3fromArray(tmpV, vertices, i * stride), t);
                StringBuilder.writeSafe(vertexBuilder, (i === 0) ? '(' : ',(');
                StringBuilder.writeFloat(vertexBuilder, tmpV[0], 10000);
                StringBuilder.writeSafe(vertexBuilder, ',');
                StringBuilder.writeFloat(vertexBuilder, tmpV[1], 10000);
                StringBuilder.writeSafe(vertexBuilder, ',');
                StringBuilder.writeFloat(vertexBuilder, tmpV[2], 10000);
                StringBuilder.writeSafe(vertexBuilder, ')');
            }

            // normal
            for (let i = 0; i < vertexCount; ++i) {
                v3transformMat3(tmpV, v3fromArray(tmpV, normals!, i * stride), n);
                StringBuilder.writeSafe(normalBuilder, (i === 0) ? '(' : ',(');
                StringBuilder.writeFloat(normalBuilder, tmpV[0], 100);
                StringBuilder.writeSafe(normalBuilder, ',');
                StringBuilder.writeFloat(normalBuilder, tmpV[1], 100);
                StringBuilder.writeSafe(normalBuilder, ',');
                StringBuilder.writeFloat(normalBuilder, tmpV[2], 100);
                StringBuilder.writeSafe(normalBuilder, ')');
            }

            const geoData = { values, groups, vertexCount, instanceIndex, isGeoTexture, mode };

            // face
            for (let i = 0; i < drawCount; ++i) {
                const v = isGeoTexture ? i : indices![i];
                if (i > 0) StringBuilder.writeSafe(indexBuilder, ',');
                StringBuilder.writeInteger(indexBuilder, v);
            }

            // color
            const quantizedColors = new Uint8Array(drawCount * 3);
            for (let i = 0; i < drawCount; i += 3) {
                const v = isGeoTexture ? i : indices![i];
                const color = UsdzExporter.getColor(v, geoData, interpolatedColors, interpolatedOverpaint);
                Color.toArray(color, quantizedColors, i);
            }
            UsdzExporter.quantizeColors(quantizedColors, vertexCount);

            // material
            const faceIndicesByMaterial = new Map<number, number[]>();
            for (let i = 0; i < drawCount; i += 3) {
                const color = Color.fromArray(quantizedColors, i);

                const transparency = UsdzExporter.getTransparency(i, geoData, interpolatedTransparency);
                const alpha = Math.round(uAlpha * (1 - transparency) * 10) / 10; // quantized

                const materialKey = this.addMaterial(color, alpha, metalness, roughness);
                let faceIndices = faceIndicesByMaterial.get(materialKey);
                if (faceIndices === undefined) {
                    faceIndices = [];
                    faceIndicesByMaterial.set(materialKey, faceIndices);
                }
                faceIndices.push(i / 3);
            }

            // If this mesh uses only one material, bind it to the material directly.
            // Otherwise, use GeomSubsets to bind it to multiple materials.
            let materialBinding: string;
            if (faceIndicesByMaterial.size === 1) {
                const materialKey = faceIndicesByMaterial.keys().next().value;
                materialBinding = `rel material:binding = </material${materialKey}>`;
            } else {
                const geomSubsets: string[] = [];
                faceIndicesByMaterial.forEach((faceIndices: number[], materialKey: number) => {
                    geomSubsets.push(`
    def GeomSubset "g${materialKey}"
    {
        uniform token elementType = "face"
        uniform token familyName = "materialBind"
        int[] indices = [${faceIndices.join(',')}]
        rel material:binding = </material${materialKey}>
    }
`);
                });
                materialBinding = geomSubsets.join('');
            }

            this.meshes.push(`
def Mesh "mesh${this.meshes.length}"
{
    int[] faceVertexCounts = [${new Array(drawCount / 3).fill(3).join(',')}]
    int[] faceVertexIndices = [${StringBuilder.getString(indexBuilder)}]
    point3f[] points = [${StringBuilder.getString(vertexBuilder)}]
    normal3f[] primvars:normals = [${StringBuilder.getString(normalBuilder)}] (
        interpolation = "vertex"
    )
    uniform token subdivisionScheme = "none"
    ${materialBinding}
}
`);
        }
    }

    async getData(ctx: RuntimeContext) {
        const header = `#usda 1.0
(
    customLayerData = {
        string creator = "Mol* ${PLUGIN_VERSION}"
    }
    metersPerUnit = 1
)
`;
        const usda = [header, ...this.materials, ...this.meshes].join('');
        const usdaData = new Uint8Array(usda.length);
        asciiWrite(usdaData, usda);
        const zipDataObj = {
            ['model.usda']: usdaData
        };
        return {
            usdz: await zip(ctx, zipDataObj, true)
        };
    }

    async getBlob(ctx: RuntimeContext) {
        const { usdz } = await this.getData(ctx);
        return new Blob([usdz], { type: 'model/vnd.usdz+zip' });
    }

    constructor(boundingBox: Box3D, radius: number) {
        super();
        const t = Mat4();
        // scale the model so that it fits within 1 meter
        Mat4.fromUniformScaling(t, Math.min(1 / (radius * 2), 1));
        // translate the model so that it sits on the ground plane (y = 0)
        Mat4.translate(t, t, Vec3.create(
            -(boundingBox.min[0] + boundingBox.max[0]) / 2,
            -boundingBox.min[1],
            -(boundingBox.min[2] + boundingBox.max[2]) / 2
        ));
        this.centerTransform = t;
    }
}