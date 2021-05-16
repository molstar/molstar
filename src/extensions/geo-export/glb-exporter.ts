/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import { asciiWrite } from '../../mol-io/common/ascii';
import { IsNativeEndianLittle, flipByteOrder } from '../../mol-io/common/binary';
import { Vec3, Mat3, Mat4 } from '../../mol-math/linear-algebra';
import { RuntimeContext } from '../../mol-task';
import { Color } from '../../mol-util/color/color';
import { arrayMinMax, fillSerial } from '../../mol-util/array';
import { NumberArray } from '../../mol-util/type-helpers';
import { MeshExporter, AddMeshInput } from './mesh-exporter';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3transformMat4 = Vec3.transformMat4;
const v3transformMat3 = Vec3.transformMat3;
const v3normalize = Vec3.normalize;
const v3toArray = Vec3.toArray;
const mat3directionTransform = Mat3.directionTransform;

// https://github.com/KhronosGroup/glTF/tree/master/specification/2.0

export type GlbData = {
    glb: Uint8Array
}

export class GlbExporter extends MeshExporter<GlbData> {
    readonly fileExtension = 'glb';
    private primitives: Record<string, any>[] = [];
    private accessors: Record<string, any>[] = [];
    private bufferViews: Record<string, any>[] = [];
    private binaryBuffer: ArrayBuffer[] = [];
    private byteOffset = 0;

    private static vecMinMax(a: NumberArray, length: number) {
        const min: number[] = (new Array(length)).fill(Infinity);
        const max: number[] = (new Array(length)).fill(-Infinity);
        for (let i = 0, il = a.length; i < il; i += length) {
            for (let j = 0; j < length; ++j) {
                min[j] = Math.min(a[i + j], min[j]);
                max[j] = Math.max(a[i + j], max[j]);
            }
        }
        return [ min, max ];
    }

    protected async addMeshWithColors(input: AddMeshInput) {
        const { mesh, meshes, values, isGeoTexture, webgl, ctx } = input;

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
            interpolatedColors = GlbExporter.getInterpolatedColors(mesh!.vertices, mesh!.vertexCount, values, stride, colorType, webgl!);
        }

        await ctx.update({ isIndeterminate: false, current: 0, max: instanceCount });

        for (let instanceIndex = 0; instanceIndex < instanceCount; ++instanceIndex) {
            if (ctx.shouldUpdate) await ctx.update({ current: instanceIndex + 1 });

            let vertices: Float32Array;
            let normals: Float32Array;
            let indices: Uint32Array | undefined;
            let groups: Float32Array | Uint8Array;
            let vertexCount: number;
            let drawCount: number;
            if (mesh !== undefined) {
                vertices = mesh.vertices;
                normals = mesh.normals;
                indices = mesh.indices;
                groups = mesh.groups;
                vertexCount = mesh.vertexCount;
                drawCount = mesh.drawCount;
            } else {
                const mesh = meshes![instanceIndex];
                vertices = mesh.vertexBuffer.ref.value;
                normals = mesh.normalBuffer.ref.value;
                indices = mesh.indexBuffer.ref.value;
                groups = mesh.groupBuffer.ref.value;
                vertexCount = mesh.vertexCount;
                drawCount = indices.length;
            }

            Mat4.fromArray(t, aTransform, instanceIndex * 16);
            mat3directionTransform(n, t);

            const vertexArray = new Float32Array(vertexCount * 3);
            const normalArray = new Float32Array(vertexCount * 3);
            const colorArray = new Float32Array(vertexCount * 4);
            let indexArray: Uint32Array;

            // position
            for (let i = 0; i < vertexCount; ++i) {
                v3transformMat4(tmpV, v3fromArray(tmpV, vertices, i * stride), t);
                v3toArray(tmpV, vertexArray, i * 3);
            }

            // normal
            for (let i = 0; i < vertexCount; ++i) {
                v3fromArray(tmpV, normals, i * stride);
                v3transformMat3(tmpV, v3normalize(tmpV, tmpV), n);
                v3toArray(tmpV, normalArray, i * 3);
            }

            // color
            for (let i = 0; i < vertexCount; ++i) {
                let color: Color;
                switch (colorType) {
                    case 'uniform':
                        color = Color.fromNormalizedArray(values.uColor.ref.value, 0);
                        break;
                    case 'instance':
                        color = Color.fromArray(tColor, instanceIndex * 3);
                        break;
                    case 'group': {
                        const group = isGeoTexture ? GlbExporter.getGroup(groups, i) : groups[i];
                        color = Color.fromArray(tColor, group * 3);
                        break;
                    }
                    case 'groupInstance': {
                        const group = isGeoTexture ? GlbExporter.getGroup(groups, i) : groups[i];
                        color = Color.fromArray(tColor, (instanceIndex * groupCount + group) * 3);
                        break;
                    }
                    case 'vertex':
                        color = Color.fromArray(tColor, i * 3);
                        break;
                    case 'vertexInstance':
                        color = Color.fromArray(tColor, (instanceIndex * drawCount + i) * 3);
                        break;
                    case 'volume':
                        color = Color.fromArray(interpolatedColors!, i * 3);
                        break;
                    case 'volumeInstance':
                        color = Color.fromArray(interpolatedColors!, (instanceIndex * vertexCount + i) * 3);
                        break;
                    default: throw new Error('Unsupported color type.');
                }

                let alpha = uAlpha;
                if (dTransparency) {
                    const group = isGeoTexture ? GlbExporter.getGroup(groups, i) : groups[i];
                    const transparency = tTransparency.array[instanceIndex * groupCount + group] / 255;
                    alpha *= 1 - transparency;
                }

                Color.toArrayNormalized(color, colorArray, i * 4);
                colorArray[i * 4 + 3] = alpha;
            }

            // face
            if (isGeoTexture) {
                indexArray = new Uint32Array(drawCount);
                fillSerial(indexArray);
            } else {
                indexArray = indices!.slice(0, drawCount);
            }

            const [ vertexMin, vertexMax ] = GlbExporter.vecMinMax(vertexArray, 3);
            const [ normalMin, normalMax ] = GlbExporter.vecMinMax(normalArray, 3);
            const [ colorMin, colorMax ] = GlbExporter.vecMinMax(colorArray, 4);
            const [ indexMin, indexMax ] = arrayMinMax(indexArray);

            // binary buffer
            let vertexBuffer = vertexArray.buffer;
            let normalBuffer = normalArray.buffer;
            let colorBuffer = colorArray.buffer;
            let indexBuffer = indexArray.buffer;
            if (!IsNativeEndianLittle) {
                vertexBuffer = flipByteOrder(new Uint8Array(vertexBuffer), 4);
                normalBuffer = flipByteOrder(new Uint8Array(normalBuffer), 4);
                colorBuffer = flipByteOrder(new Uint8Array(colorBuffer), 4);
                indexBuffer = flipByteOrder(new Uint8Array(indexBuffer), 4);
            }
            this.binaryBuffer.push(vertexBuffer, normalBuffer, colorBuffer, indexBuffer);

            // buffer views
            const bufferViewOffset = this.bufferViews.length;

            this.bufferViews.push({
                buffer: 0,
                byteOffset: this.byteOffset,
                byteLength: vertexBuffer.byteLength,
                target: 34962 // ARRAY_BUFFER
            });
            this.byteOffset += vertexBuffer.byteLength;

            this.bufferViews.push({
                buffer: 0,
                byteOffset: this.byteOffset,
                byteLength: normalBuffer.byteLength,
                target: 34962 // ARRAY_BUFFER
            });
            this.byteOffset += normalBuffer.byteLength;

            this.bufferViews.push({
                buffer: 0,
                byteOffset: this.byteOffset,
                byteLength: colorBuffer.byteLength,
                target: 34962 // ARRAY_BUFFER
            });
            this.byteOffset += colorBuffer.byteLength;

            this.bufferViews.push({
                buffer: 0,
                byteOffset: this.byteOffset,
                byteLength: indexBuffer.byteLength,
                target: 34963 // ELEMENT_ARRAY_BUFFER
            });
            this.byteOffset += indexBuffer.byteLength;

            // accessors
            const accessorOffset = this.accessors.length;
            this.accessors.push({
                bufferView: bufferViewOffset,
                byteOffset: 0,
                componentType: 5126, // FLOAT
                count: vertexCount,
                type: 'VEC3',
                max: vertexMax,
                min: vertexMin
            });
            this.accessors.push({
                bufferView: bufferViewOffset + 1,
                byteOffset: 0,
                componentType: 5126, // FLOAT
                count: vertexCount,
                type: 'VEC3',
                max: normalMax,
                min: normalMin
            });
            this.accessors.push({
                bufferView: bufferViewOffset + 2,
                byteOffset: 0,
                componentType: 5126, // FLOAT
                count: vertexCount,
                type: 'VEC4',
                max: colorMax,
                min: colorMin
            });
            this.accessors.push({
                bufferView: bufferViewOffset + 3,
                byteOffset: 0,
                componentType: 5125, // UNSIGNED_INT
                count: drawCount,
                type: 'SCALAR',
                max: [ indexMax ],
                min: [ indexMin ]
            });

            // primitive
            this.primitives.push({
                attributes: {
                    POSITION: accessorOffset,
                    NORMAL: accessorOffset + 1,
                    COLOR_0: accessorOffset + 2,
                },
                indices: accessorOffset + 3,
                material: 0
            });
        }
    }

    getData() {
        const binaryBufferLength = this.byteOffset;

        const gltf = {
            asset: {
                version: '2.0'
            },
            scenes: [{
                nodes: [ 0 ]
            }],
            nodes: [{
                mesh: 0
            }],
            meshes: [{
                primitives: this.primitives
            }],
            buffers: [{
                byteLength: binaryBufferLength,
            }],
            bufferViews: this.bufferViews,
            accessors: this.accessors,
            materials: [{}]
        };

        const createChunk = (chunkType: number, data: ArrayBuffer[], byteLength: number, padChar: number): [ArrayBuffer[], number] => {
            let padding = null;
            if (byteLength % 4 !== 0) {
                const pad = 4 - (byteLength % 4);
                byteLength += pad;
                padding = new Uint8Array(pad);
                padding.fill(padChar);
            }
            const preamble = new ArrayBuffer(8);
            const preambleDataView = new DataView(preamble);
            preambleDataView.setUint32(0, byteLength, true);
            preambleDataView.setUint32(4, chunkType, true);
            const chunk = [preamble, ...data];
            if (padding) {
                chunk.push(padding.buffer);
            }
            return [ chunk, 8 + byteLength ];
        };
        const jsonString = JSON.stringify(gltf);
        const jsonBuffer = new Uint8Array(jsonString.length);
        asciiWrite(jsonBuffer, jsonString);

        const [ jsonChunk, jsonChunkLength ] = createChunk(0x4E4F534A, [jsonBuffer.buffer], jsonBuffer.length, 0x20);
        const [ binaryChunk, binaryChunkLength ] = createChunk(0x004E4942, this.binaryBuffer, binaryBufferLength, 0x00);

        const glbBufferLength = 12 + jsonChunkLength + binaryChunkLength;
        const header = new ArrayBuffer(12);
        const headerDataView = new DataView(header);
        headerDataView.setUint32(0, 0x46546C67, true); // magic number "glTF"
        headerDataView.setUint32(4, 2, true); // version
        headerDataView.setUint32(8, glbBufferLength, true); // length
        const glbBuffer = [header, ...jsonChunk, ...binaryChunk];

        const glb = new Uint8Array(glbBufferLength);
        let offset = 0;
        for (const buffer of glbBuffer) {
            glb.set(new Uint8Array(buffer), offset);
            offset += buffer.byteLength;
        }
        return { glb };
    }

    async getBlob(ctx: RuntimeContext) {
        return new Blob([this.getData().glb], { type: 'model/gltf-binary' });
    }
}