/**
 * Copyright (c) 2021-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { asciiWrite } from '../../mol-io/common/ascii';
import { IsNativeEndianLittle, flipByteOrder } from '../../mol-io/common/binary';
import { Box3D } from '../../mol-math/geometry';
import { Vec3, Mat4 } from '../../mol-math/linear-algebra';
import { PLUGIN_VERSION } from '../../mol-plugin/version';
import { RuntimeContext } from '../../mol-task';
import { Color } from '../../mol-util/color/color';
import { fillSerial } from '../../mol-util/array';
import { NumberArray, assertUnreachable } from '../../mol-util/type-helpers';
import { MeshExporter, AddMeshInput, MeshGeoData } from './mesh-exporter';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3normalize = Vec3.normalize;
const v3toArray = Vec3.toArray;

// https://github.com/KhronosGroup/glTF/tree/master/specification/2.0

const UNSIGNED_BYTE = 5121;
const UNSIGNED_INT = 5125;
const FLOAT = 5126;
const ARRAY_BUFFER = 34962;
const ELEMENT_ARRAY_BUFFER = 34963;

const GLTF_MAGIC_BYTE = 0x46546C67;
const JSON_CHUNK_TYPE = 0x4E4F534A;
const BIN_CHUNK_TYPE = 0x004E4942;
const JSON_PAD_CHAR = 0x20;
const BIN_PAD_CHAR = 0x00;

function getPrimitiveMode(mode: 'points' | 'lines' | 'triangles'): number {
    switch (mode) {
        case 'points': return 0;
        case 'lines': return 1;
        case 'triangles': return 4;
        default: assertUnreachable(mode);
    }
}

export type GlbData = {
    glb: Uint8Array
}

export class GlbExporter extends MeshExporter<GlbData> {
    readonly fileExtension = 'glb';
    private nodes: Record<string, any>[] = [];
    private meshes: Record<string, any>[] = [];
    private materials: Record<string, any>[] = [];
    private materialMap = new Map<string, number>();
    private accessors: Record<string, any>[] = [];
    private bufferViews: Record<string, any>[] = [];
    private binaryBuffer: ArrayBuffer[] = [];
    private byteOffset = 0;
    private centerTransform: Mat4;

    private static vec3MinMax(a: NumberArray) {
        const min: number[] = [Infinity, Infinity, Infinity];
        const max: number[] = [-Infinity, -Infinity, -Infinity];
        for (let i = 0, il = a.length; i < il; i += 3) {
            for (let j = 0; j < 3; ++j) {
                min[j] = Math.min(a[i + j], min[j]);
                max[j] = Math.max(a[i + j], max[j]);
            }
        }
        return [min, max];
    }

    private addBuffer(buffer: ArrayBuffer, componentType: number, type: string, count: number, target: number, min?: any, max?: any, normalized?: boolean) {
        this.binaryBuffer.push(buffer);

        const bufferViewOffset = this.bufferViews.length;
        this.bufferViews.push({
            buffer: 0,
            byteOffset: this.byteOffset,
            byteLength: buffer.byteLength,
            target
        });
        this.byteOffset += buffer.byteLength;

        const accessorOffset = this.accessors.length;
        this.accessors.push({
            bufferView: bufferViewOffset,
            byteOffset: 0,
            componentType,
            count,
            type,
            min,
            max,
            normalized
        });
        return accessorOffset;
    }

    private addGeometryBuffers(vertices: Float32Array, normals: Float32Array | undefined, indices: Uint32Array | undefined, vertexCount: number, drawCount: number, isGeoTexture: boolean) {
        const tmpV = Vec3();
        const stride = isGeoTexture ? 4 : 3;

        const vertexArray = new Float32Array(vertexCount * 3);
        let normalArray: Float32Array | undefined;
        let indexArray: Uint32Array | undefined;

        // position
        for (let i = 0; i < vertexCount; ++i) {
            v3fromArray(tmpV, vertices, i * stride);
            v3toArray(tmpV, vertexArray, i * 3);
        }

        // normal
        if (normals) {
            normalArray = new Float32Array(vertexCount * 3);
            for (let i = 0; i < vertexCount; ++i) {
                v3fromArray(tmpV, normals, i * stride);
                v3normalize(tmpV, tmpV);
                v3toArray(tmpV, normalArray, i * 3);
            }
        }

        // face
        if (!isGeoTexture && indices) {
            indexArray = indices.slice(0, drawCount);
        }

        const [vertexMin, vertexMax] = GlbExporter.vec3MinMax(vertexArray);

        let vertexBuffer = vertexArray.buffer;
        let normalBuffer = normalArray?.buffer;
        let indexBuffer = (isGeoTexture || !indexArray) ? undefined : indexArray.buffer;
        if (!IsNativeEndianLittle) {
            vertexBuffer = flipByteOrder(new Uint8Array(vertexBuffer), 4);
            if (normalBuffer) normalBuffer = flipByteOrder(new Uint8Array(normalBuffer), 4);
            if (!isGeoTexture) indexBuffer = flipByteOrder(new Uint8Array(indexBuffer!), 4);
        }

        return {
            vertexAccessorIndex: this.addBuffer(vertexBuffer, FLOAT, 'VEC3', vertexCount, ARRAY_BUFFER, vertexMin, vertexMax),
            normalAccessorIndex: normalBuffer ? this.addBuffer(normalBuffer, FLOAT, 'VEC3', vertexCount, ARRAY_BUFFER) : undefined,
            indexAccessorIndex: (isGeoTexture || !indexBuffer) ? undefined : this.addBuffer(indexBuffer, UNSIGNED_INT, 'SCALAR', drawCount, ELEMENT_ARRAY_BUFFER)
        };
    }

    private addColorBuffer(geoData: MeshGeoData, interpolatedColors: Uint8Array | undefined, interpolatedOverpaint: Uint8Array | undefined, interpolatedTransparency: Uint8Array | undefined) {
        const { values, vertexCount } = geoData;
        const uAlpha = values.uAlpha.ref.value;

        const colorArray = new Uint8Array(vertexCount * 4);

        for (let i = 0; i < vertexCount; ++i) {
            let color = GlbExporter.getColor(i, geoData, interpolatedColors, interpolatedOverpaint);

            const transparency = GlbExporter.getTransparency(i, geoData, interpolatedTransparency);
            const alpha = uAlpha * (1 - transparency);

            color = Color.sRGBToLinear(color);
            Color.toArray(color, colorArray, i * 4);
            colorArray[i * 4 + 3] = Math.round(alpha * 255);
        }

        let colorBuffer = colorArray.buffer;
        if (!IsNativeEndianLittle) {
            colorBuffer = flipByteOrder(new Uint8Array(colorBuffer), 4);
        }

        return this.addBuffer(colorBuffer, UNSIGNED_BYTE, 'VEC4', vertexCount, ARRAY_BUFFER, undefined, undefined, true);
    }

    private addMaterial(metalness: number, roughness: number, emissive: number, doubleSided: boolean, alpha: boolean) {
        const hash = `${metalness}|${roughness}|${emissive}|${doubleSided}`;
        if (!this.materialMap.has(hash)) {
            this.materialMap.set(hash, this.materials.length);
            this.materials.push({
                pbrMetallicRoughness: {
                    baseColorFactor: [1, 1, 1, 1],
                    metallicFactor: metalness,
                    roughnessFactor: roughness
                },
                emissiveFactor: [emissive, emissive, emissive],
                doubleSided,
                alphaMode: alpha ? 'BLEND' : 'OPAQUE',
            });
        }
        return this.materialMap.get(hash)!;
    }

    protected async addMeshWithColors(input: AddMeshInput) {
        const { mesh, values, isGeoTexture, mode, webgl, ctx } = input;

        const t = Mat4();

        const colorType = values.dColorType.ref.value;
        const overpaintType = values.dOverpaintType.ref.value;
        const transparencyType = values.dTransparencyType.ref.value;
        const dTransparency = values.dTransparency.ref.value;
        const aTransform = values.aTransform.ref.value;
        const instanceCount = values.uInstanceCount.ref.value;
        const metalness = values.uMetalness.ref.value;
        const roughness = values.uRoughness.ref.value;
        const emissive = values.uEmissive.ref.value;
        const doubleSided = values.uDoubleSided?.ref.value || values.hasReflection.ref.value;
        const alpha = values.uAlpha.ref.value < 1;

        const material = this.addMaterial(metalness, roughness, emissive, doubleSided, alpha);

        let interpolatedColors: Uint8Array | undefined;
        if (webgl && mesh && (colorType === 'volume' || colorType === 'volumeInstance')) {
            const stride = isGeoTexture ? 4 : 3;
            interpolatedColors = GlbExporter.getInterpolatedColors(webgl, { vertices: mesh.vertices, vertexCount: mesh.vertexCount, values, stride, colorType });
        }

        let interpolatedOverpaint: Uint8Array | undefined;
        if (webgl && mesh && overpaintType === 'volumeInstance') {
            const stride = isGeoTexture ? 4 : 3;
            interpolatedOverpaint = GlbExporter.getInterpolatedOverpaint(webgl, { vertices: mesh.vertices, vertexCount: mesh.vertexCount, values, stride, colorType: overpaintType });
        }

        let interpolatedTransparency: Uint8Array | undefined;
        if (webgl && mesh && transparencyType === 'volumeInstance') {
            const stride = isGeoTexture ? 4 : 3;
            interpolatedTransparency = GlbExporter.getInterpolatedTransparency(webgl, { vertices: mesh.vertices, vertexCount: mesh.vertexCount, values, stride, colorType: transparencyType });
        }

        // instancing
        const sameGeometryBuffers = mesh !== undefined;
        const sameColorBuffer = sameGeometryBuffers && colorType !== 'instance' && !colorType.endsWith('Instance') && !dTransparency;
        let vertexAccessorIndex: number;
        let normalAccessorIndex: number | undefined;
        let indexAccessorIndex: number | undefined;
        let colorAccessorIndex: number;
        let meshIndex: number;

        await ctx.update({ isIndeterminate: false, current: 0, max: instanceCount });

        for (let instanceIndex = 0; instanceIndex < instanceCount; ++instanceIndex) {
            if (ctx.shouldUpdate) await ctx.update({ current: instanceIndex + 1 });

            // create a glTF mesh if needed
            if (instanceIndex === 0 || !sameGeometryBuffers || !sameColorBuffer) {
                const { vertices, normals, indices, groups, vertexCount, drawCount } = GlbExporter.getInstance(input, instanceIndex);

                // create geometry buffers if needed
                if (instanceIndex === 0 || !sameGeometryBuffers) {
                    const accessorIndices = this.addGeometryBuffers(vertices, normals, indices, vertexCount, drawCount, isGeoTexture);
                    vertexAccessorIndex = accessorIndices.vertexAccessorIndex;
                    normalAccessorIndex = accessorIndices.normalAccessorIndex;
                    indexAccessorIndex = accessorIndices.indexAccessorIndex;
                }

                // create a color buffer if needed
                if (instanceIndex === 0 || !sameColorBuffer) {
                    colorAccessorIndex = this.addColorBuffer({ values, groups, vertexCount, instanceIndex, isGeoTexture, mode }, interpolatedColors, interpolatedOverpaint, interpolatedTransparency);
                }

                // glTF mesh
                meshIndex = this.meshes.length;
                this.meshes.push({
                    primitives: [{
                        attributes: {
                            POSITION: vertexAccessorIndex!,
                            NORMAL: normalAccessorIndex!,
                            COLOR_0: colorAccessorIndex!
                        },
                        indices: indexAccessorIndex,
                        material,
                        mode: getPrimitiveMode(mode),
                    }]
                });
            }

            // node
            Mat4.fromArray(t, aTransform, instanceIndex * 16);
            Mat4.mul(t, this.centerTransform, t);
            const node: Record<string, any> = {
                mesh: meshIndex!,
                matrix: t.slice()
            };
            this.nodes.push(node);
        }
    }

    async getData() {
        const binaryBufferLength = this.byteOffset;

        const gltf = {
            asset: {
                version: '2.0',
                generator: `Mol* ${PLUGIN_VERSION}`
            },
            scenes: [{
                nodes: fillSerial(new Array(this.nodes.length) as number[])
            }],
            nodes: this.nodes,
            meshes: this.meshes,
            buffers: [{
                byteLength: binaryBufferLength,
            }],
            bufferViews: this.bufferViews,
            accessors: this.accessors,
            materials: this.materials
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
            return [chunk, 8 + byteLength];
        };
        const jsonString = JSON.stringify(gltf);
        const jsonBuffer = new Uint8Array(jsonString.length);
        asciiWrite(jsonBuffer, jsonString);

        const [jsonChunk, jsonChunkLength] = createChunk(JSON_CHUNK_TYPE, [jsonBuffer.buffer], jsonBuffer.length, JSON_PAD_CHAR);
        const [binaryChunk, binaryChunkLength] = createChunk(BIN_CHUNK_TYPE, this.binaryBuffer, binaryBufferLength, BIN_PAD_CHAR);

        const glbBufferLength = 12 + jsonChunkLength + binaryChunkLength;
        const header = new ArrayBuffer(12);
        const headerDataView = new DataView(header);
        headerDataView.setUint32(0, GLTF_MAGIC_BYTE, true); // magic number "glTF"
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
        return new Blob([(await this.getData()).glb], { type: 'model/gltf-binary' });
    }

    constructor(boundingBox: Box3D) {
        super();
        const tmpV = Vec3();
        Vec3.add(tmpV, boundingBox.min, boundingBox.max);
        Vec3.scale(tmpV, tmpV, -0.5);
        this.centerTransform = Mat4.fromTranslation(Mat4(), tmpV);
    }
}