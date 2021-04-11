/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { MeshValues } from '../../mol-gl/renderable/mesh';
import { LinesValues } from '../../mol-gl/renderable/lines';
import { PointsValues } from '../../mol-gl/renderable/points';
import { SpheresValues } from '../../mol-gl/renderable/spheres';
import { CylindersValues } from '../../mol-gl/renderable/cylinders';
import { BaseValues, SizeValues } from '../../mol-gl/renderable/schema';
import { TextureImage } from '../../mol-gl/renderable/util';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
import { addCylinder } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { Vec3, Mat3, Mat4 } from '../../mol-math/linear-algebra';
import { RuntimeContext } from '../../mol-task';
import { StringBuilder } from '../../mol-util';
import { Color } from '../../mol-util/color/color';
import { decodeFloatRGB } from '../../mol-util/float-packing';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3transformMat4 = Vec3.transformMat4;
const v3transformMat3 = Vec3.transformMat3;
const mat3directionTransform = Mat3.directionTransform;

type RenderObjectExportData = {
    [k: string]: string | Uint8Array | undefined
}

interface RenderObjectExporter<D extends RenderObjectExportData> {
    add(renderObject: GraphicsRenderObject, ctx: RuntimeContext): Promise<void>
    getData(): D
}

// http://paulbourke.net/dataformats/obj/
// http://paulbourke.net/dataformats/mtl/

export type ObjData = {
    obj: string
    mtl: string
}

export class ObjExporter implements RenderObjectExporter<ObjData> {
    private obj = StringBuilder.create();
    private mtl = StringBuilder.create();
    private vertexOffset = 0;
    private currentColor: Color | undefined;
    private currentAlpha: number | undefined;
    private materialSet = new Set<string>();

    private static getSizeFromTexture(tSize: TextureImage<Uint8Array>, i: number): number {
        const r = tSize.array[i * 3];
        const g = tSize.array[i * 3 + 1];
        const b = tSize.array[i * 3 + 2];
        return decodeFloatRGB(r, g, b);
    }

    private static getSize(values: BaseValues & SizeValues, instanceIndex: number, group: number): number {
        const tSize = values.tSize.ref.value;
        let size = 0;
        switch (values.dSizeType.ref.value) {
            case 'uniform':
                size = values.uSize.ref.value;
                break;
            case 'instance':
                size = ObjExporter.getSizeFromTexture(tSize, instanceIndex) / 100;
                break;
            case 'group':
                size = ObjExporter.getSizeFromTexture(tSize, group) / 100;
                break;
            case 'groupInstance':
                const groupCount = values.uGroupCount.ref.value;
                size = ObjExporter.getSizeFromTexture(tSize, instanceIndex * groupCount + group) / 100;
                break;
        }
        return size * values.uSizeFactor.ref.value;
    }

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

    private async addMeshWithColors(vertices: Float32Array, normals: Float32Array, indices: Uint32Array, groups: Float32Array, vertexCount: number, drawCount: number, values: BaseValues, instanceIndex: number, ctx: RuntimeContext) {
        const obj = this.obj;
        const t = Mat4();
        const n = Mat3();
        const tmpV = Vec3();

        const colorType = values.dColorType.ref.value;
        const tColor = values.tColor.ref.value.array;
        const uAlpha = values.uAlpha.ref.value;
        const aTransform = values.aTransform.ref.value;

        Mat4.fromArray(t, aTransform, instanceIndex * 16);
        mat3directionTransform(n, t);

        // position
        for (let i = 0; i < vertexCount; ++i) {
            if (i % 1000 === 0 && ctx.shouldUpdate) await ctx.update();
            v3transformMat4(tmpV, v3fromArray(tmpV, vertices, i * 3), t);
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
            if (i % 1000 === 0 && ctx.shouldUpdate) await ctx.update();
            v3transformMat3(tmpV, v3fromArray(tmpV, normals, i * 3), n);
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
            if (i % 3000 === 0 && ctx.shouldUpdate) await ctx.update();
            let color: Color;
            switch (colorType) {
                case 'uniform':
                    color = Color.fromNormalizedArray(values.uColor.ref.value, 0);
                    break;
                case 'instance':
                    color = Color.fromArray(tColor, instanceIndex * 3);
                    break;
                case 'group':
                    color = Color.fromArray(tColor, groups[indices[i]] * 3);
                    break;
                case 'groupInstance':
                    const groupCount = values.uGroupCount.ref.value;
                    const group = groups[indices[i]];
                    color = Color.fromArray(tColor, (instanceIndex * groupCount + group) * 3);
                    break;
                case 'vertex':
                    color = Color.fromArray(tColor, i * 3);
                    break;
                case 'vertexInstance':
                    color = Color.fromArray(tColor, (instanceIndex * drawCount + i) * 3);
                    break;
                default: throw new Error('Unsupported color type.');
            }
            this.updateMaterial(color, uAlpha);

            const v1 = this.vertexOffset + indices[i] + 1;
            const v2 = this.vertexOffset + indices[i + 1] + 1;
            const v3 = this.vertexOffset + indices[i + 2] + 1;
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

    private async addMesh(values: MeshValues, ctx: RuntimeContext) {
        const aPosition = values.aPosition.ref.value;
        const aNormal = values.aNormal.ref.value;
        const elements = values.elements.ref.value;
        const aGroup = values.aGroup.ref.value;
        const instanceCount = values.instanceCount.ref.value;
        const vertexCount = values.uVertexCount.ref.value;
        const drawCount = values.drawCount.ref.value;

        for (let instanceIndex = 0; instanceIndex < instanceCount; ++instanceIndex) {
            await this.addMeshWithColors(aPosition, aNormal, elements, aGroup, vertexCount, drawCount, values, instanceIndex, ctx);
        }
    }

    private async addLines(values: LinesValues, ctx: RuntimeContext) {
        // TODO
    }

    private async addPoints(values: PointsValues, ctx: RuntimeContext) {
        // TODO
    }

    private async addSpheres(values: SpheresValues, ctx: RuntimeContext) {
        const center = Vec3();

        const aPosition = values.aPosition.ref.value;
        const aGroup = values.aGroup.ref.value;
        const instanceCount = values.instanceCount.ref.value;
        const vertexCount = values.uVertexCount.ref.value;

        for (let instanceIndex = 0; instanceIndex < instanceCount; ++instanceIndex) {
            const state = MeshBuilder.createState(512, 256);

            for (let i = 0; i < vertexCount; i += 4) {
                v3fromArray(center, aPosition, i * 3);

                const group = aGroup[i];
                const radius = ObjExporter.getSize(values, instanceIndex, group);
                state.currentGroup = group;
                addSphere(state, center, radius, 2);
            }

            const mesh = MeshBuilder.getMesh(state);
            const vertices = mesh.vertexBuffer.ref.value;
            const normals = mesh.normalBuffer.ref.value;
            const indices = mesh.indexBuffer.ref.value;
            const groups = mesh.groupBuffer.ref.value;
            await this.addMeshWithColors(vertices, normals, indices, groups, vertices.length / 3, indices.length, values, instanceIndex, ctx);
        }
    }

    private async addCylinders(values: CylindersValues, ctx: RuntimeContext) {
        const start = Vec3();
        const end = Vec3();

        const aStart = values.aStart.ref.value;
        const aEnd = values.aEnd.ref.value;
        const aScale = values.aScale.ref.value;
        const aCap = values.aCap.ref.value;
        const aGroup = values.aGroup.ref.value;
        const instanceCount = values.instanceCount.ref.value;
        const vertexCount = values.uVertexCount.ref.value;

        for (let instanceIndex = 0; instanceIndex < instanceCount; ++instanceIndex) {
            const state = MeshBuilder.createState(512, 256);

            for (let i = 0; i < vertexCount; i += 6) {
                v3fromArray(start, aStart, i * 3);
                v3fromArray(end, aEnd, i * 3);

                const group = aGroup[i];
                const radius = ObjExporter.getSize(values, instanceIndex, group) * aScale[i];
                const cap = aCap[i];
                const topCap = cap === 1 || cap === 3;
                const bottomCap = cap >= 2;
                const cylinderProps = { radiusTop: radius, radiusBottom: radius, topCap, bottomCap, radialSegments: 32 };
                state.currentGroup = aGroup[i];
                addCylinder(state, start, end, 1, cylinderProps);
            }

            const mesh = MeshBuilder.getMesh(state);
            const vertices = mesh.vertexBuffer.ref.value;
            const normals = mesh.normalBuffer.ref.value;
            const indices = mesh.indexBuffer.ref.value;
            const groups = mesh.groupBuffer.ref.value;
            await this.addMeshWithColors(vertices, normals, indices, groups, vertices.length / 3, indices.length, values, instanceIndex, ctx);
        }
    }

    async add(renderObject: GraphicsRenderObject, ctx: RuntimeContext) {
        if (!renderObject.state.visible) return;

        switch (renderObject.type) {
            case 'mesh':
                await this.addMesh(renderObject.values as MeshValues, ctx);
                break;
            case 'lines':
                await this.addLines(renderObject.values as LinesValues, ctx);
                break;
            case 'points':
                await this.addPoints(renderObject.values as PointsValues, ctx);
                break;
            case 'spheres':
                await this.addSpheres(renderObject.values as SpheresValues, ctx);
                break;
            case 'cylinders':
                await this.addCylinders(renderObject.values as CylindersValues, ctx);
                break;
        }
    }

    getData() {
        return {
            obj: StringBuilder.getString(this.obj),
            mtl: StringBuilder.getString(this.mtl)
        };
    }

    constructor(filename: string) {
        StringBuilder.writeSafe(this.obj, `mtllib ${filename}.mtl\n`);
    }
}