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
import { TextureMeshValues } from '../../mol-gl/renderable/texture-mesh';
import { BaseValues, SizeValues } from '../../mol-gl/renderable/schema';
import { TextureImage } from '../../mol-gl/renderable/util';
import { WebGLContext } from '../../mol-gl/webgl/context';
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
    add(renderObject: GraphicsRenderObject, webgl: WebGLContext, ctx: RuntimeContext): Promise<void> | undefined
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

    private static getGroup(groups: Float32Array | Uint8Array, i: number): number {
        const i4 = i * 4;
        const r = groups[i4];
        const g = groups[i4 + 1];
        const b = groups[i4 + 2];
        if (groups instanceof Float32Array) {
            return decodeFloatRGB(r * 255 + 0.5, g * 255 + 0.5, b * 255 + 0.5);
        }
        return decodeFloatRGB(r, g, b);
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

    private async addMeshWithColors(vertices: Float32Array, normals: Float32Array, indices: Uint32Array | undefined, groups: Float32Array | Uint8Array, vertexCount: number, drawCount: number, values: BaseValues, instanceIndex: number, geoTexture: boolean, ctx: RuntimeContext) {
        const obj = this.obj;
        const t = Mat4();
        const n = Mat3();
        const tmpV = Vec3();
        const stride = geoTexture ? 4 : 3;

        const colorType = values.dColorType.ref.value;
        const tColor = values.tColor.ref.value.array;
        const uAlpha = values.uAlpha.ref.value;
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
                    const group = geoTexture ? ObjExporter.getGroup(groups, i) : groups[indices![i]];
                    color = Color.fromArray(tColor, group * 3);
                    break;
                }
                case 'groupInstance': {
                    const groupCount = values.uGroupCount.ref.value;
                    const group = geoTexture ? ObjExporter.getGroup(groups, i) : groups[indices![i]];
                    color = Color.fromArray(tColor, (instanceIndex * groupCount + group) * 3);
                    break;
                }
                case 'vertex':
                    color = Color.fromArray(tColor, i * 3);
                    break;
                case 'vertexInstance':
                    color = Color.fromArray(tColor, (instanceIndex * drawCount + i) * 3);
                    break;
                default: throw new Error('Unsupported color type.');
            }
            this.updateMaterial(color, uAlpha);

            const v1 = this.vertexOffset + (geoTexture ? i : indices![i]) + 1;
            const v2 = this.vertexOffset + (geoTexture ? i + 1 : indices![i + 1]) + 1;
            const v3 = this.vertexOffset + (geoTexture ? i + 2 : indices![i + 2]) + 1;
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
            await this.addMeshWithColors(aPosition, aNormal, elements, aGroup, vertexCount, drawCount, values, instanceIndex, false, ctx);
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
            await this.addMeshWithColors(vertices, normals, indices, groups, vertices.length / 3, indices.length, values, instanceIndex, false, ctx);
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
            await this.addMeshWithColors(vertices, normals, indices, groups, vertices.length / 3, indices.length, values, instanceIndex, false, ctx);
        }
    }

    private async addTextureMesh(values: TextureMeshValues, webgl: WebGLContext, ctx: RuntimeContext) {
        const GeoExportName = 'geo-export';
        if (!webgl.namedFramebuffers[GeoExportName]) {
            webgl.namedFramebuffers[GeoExportName] = webgl.resources.framebuffer();
        }
        const framebuffer = webgl.namedFramebuffers[GeoExportName];

        const [ width, height ] = values.uGeoTexDim.ref.value;
        const vertices = new Float32Array(width * height * 4);
        const normals = new Float32Array(width * height * 4);
        const groups = webgl.isWebGL2 ? new Uint8Array(width * height * 4) : new Float32Array(width * height * 4);

        framebuffer.bind();
        values.tPosition.ref.value.attachFramebuffer(framebuffer, 0);
        webgl.readPixels(0, 0, width, height, vertices);
        values.tNormal.ref.value.attachFramebuffer(framebuffer, 0);
        webgl.readPixels(0, 0, width, height, normals);
        values.tGroup.ref.value.attachFramebuffer(framebuffer, 0);
        webgl.readPixels(0, 0, width, height, groups);

        const vertexCount = values.uVertexCount.ref.value;
        const instanceCount = values.instanceCount.ref.value;
        const drawCount = values.drawCount.ref.value;

        for (let instanceIndex = 0; instanceIndex < instanceCount; ++instanceIndex) {
            await this.addMeshWithColors(vertices, normals, undefined, groups, vertexCount, drawCount, values, instanceIndex, true, ctx);
        }
    }

    add(renderObject: GraphicsRenderObject, webgl: WebGLContext, ctx: RuntimeContext) {
        if (!renderObject.state.visible) return;

        switch (renderObject.type) {
            case 'mesh':
                return this.addMesh(renderObject.values as MeshValues, ctx);
            case 'lines':
                return this.addLines(renderObject.values as LinesValues, ctx);
            case 'points':
                return this.addPoints(renderObject.values as PointsValues, ctx);
            case 'spheres':
                return this.addSpheres(renderObject.values as SpheresValues, ctx);
            case 'cylinders':
                return this.addCylinders(renderObject.values as CylindersValues, ctx);
            case 'texture-mesh':
                return this.addTextureMesh(renderObject.values as TextureMeshValues, webgl, ctx);
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