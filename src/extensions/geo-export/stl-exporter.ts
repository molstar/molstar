/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import { BaseValues } from '../../mol-gl/renderable/schema';
import { Vec3, Mat4 } from '../../mol-math/linear-algebra';
import { RuntimeContext } from '../../mol-task';
import { MeshExporter } from './mesh-exporter';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3transformMat4 = Vec3.transformMat4;
const v3triangleNormal = Vec3.triangleNormal;
const v3toArray = Vec3.toArray;

export type StlData = {
    stl: Uint8Array
}

export class StlExporter extends MeshExporter<StlData> {
    readonly fileExtension = 'stl';
    private triangleBuffers: ArrayBuffer[] = [];
    private triangleCount = 0;

    async addMeshWithColors(vertices: Float32Array, normals: Float32Array, indices: Uint32Array | undefined, groups: Float32Array | Uint8Array, vertexCount: number, drawCount: number, values: BaseValues, instanceIndex: number, geoTexture: boolean, ctx: RuntimeContext) {
        const t = Mat4();
        const tmpV = Vec3();
        const v1 = Vec3();
        const v2 = Vec3();
        const v3 = Vec3();
        const n = Vec3();
        const stride = geoTexture ? 4 : 3;

        const aTransform = values.aTransform.ref.value;
        Mat4.fromArray(t, aTransform, instanceIndex * 16);

        const currentProgress = (vertexCount + drawCount) * instanceIndex;
        await ctx.update({ isIndeterminate: false, current: currentProgress, max: (vertexCount + drawCount) * values.uInstanceCount.ref.value });

        // position
        const vertexArray = new Float32Array(vertexCount * 3);
        for (let i = 0; i < vertexCount; ++i) {
            if (i % 1000 === 0 && ctx.shouldUpdate) await ctx.update({ current: currentProgress + i });
            v3transformMat4(tmpV, v3fromArray(tmpV, vertices, i * stride), t);
            v3toArray(tmpV, vertexArray, i * 3);
        }

        // face
        const triangleBuffer = new ArrayBuffer(50 * drawCount);
        const dataView = new DataView(triangleBuffer);
        for (let i = 0; i < drawCount; i += 3) {
            if (i % 3000 === 0 && ctx.shouldUpdate) await ctx.update({ current: currentProgress + vertexCount + i });

            v3fromArray(v1, vertexArray, (geoTexture ? i : indices![i]) * 3);
            v3fromArray(v2, vertexArray, (geoTexture ? i + 1 : indices![i + 1]) * 3);
            v3fromArray(v3, vertexArray, (geoTexture ? i + 2 : indices![i + 2]) * 3);
            v3triangleNormal(n, v1, v2, v3);

            const byteOffset = 50 * i;
            dataView.setFloat32(byteOffset, n[0], true);
            dataView.setFloat32(byteOffset + 4, n[1], true);
            dataView.setFloat32(byteOffset + 8, n[2], true);

            dataView.setFloat32(byteOffset + 12, v1[0], true);
            dataView.setFloat32(byteOffset + 16, v1[1], true);
            dataView.setFloat32(byteOffset + 20, v1[2], true);

            dataView.setFloat32(byteOffset + 24, v2[0], true);
            dataView.setFloat32(byteOffset + 28, v2[1], true);
            dataView.setFloat32(byteOffset + 32, v2[2], true);

            dataView.setFloat32(byteOffset + 36, v3[0], true);
            dataView.setFloat32(byteOffset + 40, v3[1], true);
            dataView.setFloat32(byteOffset + 44, v3[2], true);
        }

        this.triangleBuffers.push(triangleBuffer);
        this.triangleCount += drawCount;
    }

    getData() {
        const stl = new Uint8Array(84 + 50 * this.triangleCount);

        const dataView = new DataView(stl.buffer);
        dataView.setUint32(80, this.triangleCount, true);

        let byteOffset = 84;
        for (const buffer of this.triangleBuffers) {
            stl.set(new Uint8Array(buffer), byteOffset);
            byteOffset += buffer.byteLength;
        }
        return { stl };
    }

    async getBlob(ctx: RuntimeContext) {
        return new Blob([this.getData().stl], { type: 'model/stl' });
    }
}