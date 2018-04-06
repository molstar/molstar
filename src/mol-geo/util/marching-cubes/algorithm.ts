/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task, RuntimeContext } from 'mol-task'
import { ChunkedArray } from 'mol-data/util'
import { Tensor } from 'mol-math/linear-algebra'
import { Surface } from '../../shape/surface'
import { Index, EdgeIdInfo, CubeEdges, EdgeTable, TriTable } from './tables'
import { ValueCell } from 'mol-util'

/**
 * The parameters required by the algorithm.
 */
export interface MarchingCubesParameters {
    isoLevel: number,
    scalarField: Tensor,

    bottomLeft?: ArrayLike<number>,
    topRight?: ArrayLike<number>,

    annotationField?: Tensor,

    oldSurface?: Surface
}

export function compute(parameters: MarchingCubesParameters) {
    return Task.create('Marching Cubes', async ctx => {
        let comp = new MarchingCubesComputation(parameters, ctx);
        return await comp.run();
    });
}

class MarchingCubesComputation {
    private size: number;
    private sliceSize: number;
    private parameters: MarchingCubesParameters;

    private minX = 0; private minY = 0; private minZ = 0;
    private maxX = 0; private maxY = 0; private maxZ = 0;
    private state: MarchingCubesState;

    private async doSlices() {
        let done = 0;

        for (let k = this.minZ; k < this.maxZ; k++) {
            this.slice(k);

            done += this.sliceSize;
            if (this.ctx.shouldUpdate) {
                await this.ctx.update({ message: 'Computing surface...', current: done, max: this.size });
            }
        }
    }

    private slice(k: number) {
        for (let j = this.minY; j < this.maxY; j++) {
            for (let i = this.minX; i < this.maxX; i++) {
                this.state.processCell(i, j, k);
            }
        }
        this.state.clearEdgeVertexIndexSlice(k);
    }

    private finish() {
        const vb = ChunkedArray.compact(this.state.vertexBuffer, true) as Float32Array;
        const ib = ChunkedArray.compact(this.state.triangleBuffer, true) as Uint32Array;

        this.state.vertexBuffer = <any>void 0;
        this.state.verticesOnEdges = <any>void 0;

        let ret: Surface = {
            vertexCount:  this.state.vertexCount,
            triangleCount: this.state.triangleCount,
            vertexBuffer: this.parameters.oldSurface ? ValueCell.update(this.parameters.oldSurface.vertexBuffer, vb) : ValueCell.create(vb),
            indexBuffer: this.parameters.oldSurface ? ValueCell.update(this.parameters.oldSurface.indexBuffer, ib) : ValueCell.create(ib),
            normalBuffer: this.parameters.oldSurface ? this.parameters.oldSurface.normalBuffer : ValueCell.create(void 0),
            vertexAnnotation: this.state.annotate
                ? this.parameters.oldSurface && this.parameters.oldSurface.vertexAnnotation
                    ? ValueCell.update(this.parameters.oldSurface.vertexAnnotation, ChunkedArray.compact(this.state.annotationBuffer))
                    : ValueCell.create(ChunkedArray.compact(this.state.annotationBuffer))
                : void 0,
            normalsComputed: false
        }

        return ret;
    }

    async run() {
        await this.ctx.update({ message: 'Computing surface...', current: 0, max: this.size });
        await this.doSlices();
        await this.ctx.update('Finalizing...');
        return this.finish();
    }

    constructor(
        parameters: MarchingCubesParameters,
        private ctx: RuntimeContext) {

        let params = { ...parameters };
        this.parameters = params;

        if (!params.bottomLeft) params.bottomLeft = [0, 0, 0];
        if (!params.topRight) params.topRight = params.scalarField.space.dimensions;

        this.state = new MarchingCubesState(params),
            this.minX = params.bottomLeft[0]; this.minY = params.bottomLeft[1]; this.minZ = params.bottomLeft[2];
        this.maxX = params.topRight[0] - 1; this.maxY = params.topRight[1] - 1; this.maxZ = params.topRight[2] - 1;

        this.size = (this.maxX - this.minX) * (this.maxY - this.minY) * (this.maxZ - this.minZ);
        this.sliceSize = (this.maxX - this.minX) * (this.maxY - this.minY);
    }
}

class MarchingCubesState {
    nX: number; nY: number; nZ: number;
    isoLevel: number;
    scalarFieldGet: Tensor.Space['get'];
    scalarField: Tensor.Data;
    annotationFieldGet?: Tensor.Space['get'];
    annotationField?: Tensor.Data;
    annotate: boolean;

    // two layers of vertex indices. Each vertex has 3 edges associated.
    verticesOnEdges: Int32Array;
    vertList: number[] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    i: number = 0; j: number = 0; k: number = 0;

    vertexBuffer: ChunkedArray<number>;
    annotationBuffer: ChunkedArray<number>;
    triangleBuffer: ChunkedArray<number>;
    vertexCount = 0;
    triangleCount = 0;

    private get3dOffsetFromEdgeInfo(index: Index) {
        return (this.nX * (((this.k + index.k) % 2) * this.nY + this.j + index.j) + this.i + index.i);
    }

    /**
     * This clears the "vertex index buffer" for the slice that will not be accessed anymore.
     */
    clearEdgeVertexIndexSlice(k: number) {
        // clear either the top or bottom half of the buffer...
        const start = k % 2 === 0 ? 0 : 3 * this.nX * this.nY;
        const end = k % 2 === 0 ? 3 * this.nX * this.nY : this.verticesOnEdges.length;
        for (let i = start; i < end; i++) this.verticesOnEdges[i] = 0;
    }

    private interpolate(edgeNum: number) {
        const info = EdgeIdInfo[edgeNum],
            edgeId = 3 * this.get3dOffsetFromEdgeInfo(info) + info.e;

        const ret = this.verticesOnEdges[edgeId];
        if (ret > 0) return (ret - 1) | 0;

        const edge = CubeEdges[edgeNum];
        const a = edge.a, b = edge.b;
        const li = a.i + this.i, lj = a.j + this.j, lk = a.k + this.k;
        const hi = b.i + this.i, hj = b.j + this.j, hk = b.k + this.k;
        const v0 = this.scalarFieldGet(this.scalarField, li, lj, lk), v1 = this.scalarFieldGet(this.scalarField, hi, hj, hk);
        const t = (this.isoLevel - v0) / (v0 - v1);

        const id = ChunkedArray.add3(
            this.vertexBuffer,
            li + t * (li - hi),
            lj + t * (lj - hj),
            lk + t * (lk - hk)) | 0;

        this.verticesOnEdges[edgeId] = id + 1;

        if (this.annotate) {
            const u = this.annotationFieldGet!(this.annotationField!, li, lj, lk);
            const v = this.annotationFieldGet!(this.annotationField!, hi, hj, hk)
            let a = t < 0.5 ? u : v;
            if (a < 0) a = t < 0.5 ? v : u;
            ChunkedArray.add(this.annotationBuffer, a);
        }

        this.vertexCount++;

        return id;
    }

    constructor(params: MarchingCubesParameters) {
        const dims = params.scalarField.space.dimensions;
        this.nX = dims[0]; this.nY = dims[1]; this.nZ = dims[2];
        this.isoLevel = params.isoLevel;
        this.scalarFieldGet = params.scalarField.space.get;
        this.scalarField = params.scalarField.data;
        if (params.annotationField) {
            this.annotationField = params.annotationField.data;
            this.annotationFieldGet = params.annotationField.space.get;
        }

        let dX = params.topRight![0] - params.bottomLeft![0], dY = params.topRight![1] - params.bottomLeft![1], dZ = params.topRight![2] - params.bottomLeft![2],
            vertexBufferSize = Math.min(262144, Math.max(dX * dY * dZ / 16, 1024) | 0),
            triangleBufferSize = Math.min(1 << 16, vertexBufferSize * 4);

        this.vertexBuffer = ChunkedArray.create<number>(s => new Float32Array(s), 3, vertexBufferSize,
            params.oldSurface && params.oldSurface.vertexBuffer.ref.value);
        this.triangleBuffer = ChunkedArray.create<number>(s => new Uint32Array(s), 3, triangleBufferSize,
            params.oldSurface && params.oldSurface.indexBuffer.ref.value);

        this.annotate = !!params.annotationField;
        if (this.annotate) this.annotationBuffer = ChunkedArray.create(s => new Int32Array(s), 1, vertexBufferSize);

        // two layers of vertex indices. Each vertex has 3 edges associated.
        this.verticesOnEdges = new Int32Array(3 * this.nX * this.nY * 2);
    }

    get(i: number, j: number, k: number) {
        return this.scalarFieldGet(this.scalarField, i, j, k);
    }

    processCell(i: number, j: number, k: number) {
        let tableIndex = 0;

        if (this.get(i, j, k) < this.isoLevel) tableIndex |= 1;
        if (this.get(i + 1, j, k) < this.isoLevel) tableIndex |= 2;
        if (this.get(i + 1, j + 1, k) < this.isoLevel) tableIndex |= 4;
        if (this.get(i, j + 1, k) < this.isoLevel) tableIndex |= 8;
        if (this.get(i, j, k + 1) < this.isoLevel) tableIndex |= 16;
        if (this.get(i + 1, j, k + 1) < this.isoLevel) tableIndex |= 32;
        if (this.get(i + 1, j + 1, k + 1) < this.isoLevel) tableIndex |= 64;
        if (this.get(i, j + 1, k + 1) < this.isoLevel) tableIndex |= 128;

        if (tableIndex === 0 || tableIndex === 255) return;

        this.i = i; this.j = j; this.k = k;
        let edgeInfo = EdgeTable[tableIndex];
        if ((edgeInfo & 1) > 0) this.vertList[0] = this.interpolate(0); // 0 1
        if ((edgeInfo & 2) > 0) this.vertList[1] = this.interpolate(1); // 1 2
        if ((edgeInfo & 4) > 0) this.vertList[2] = this.interpolate(2); // 2 3
        if ((edgeInfo & 8) > 0) this.vertList[3] = this.interpolate(3); // 0 3
        if ((edgeInfo & 16) > 0) this.vertList[4] = this.interpolate(4); // 4 5
        if ((edgeInfo & 32) > 0) this.vertList[5] = this.interpolate(5); // 5 6
        if ((edgeInfo & 64) > 0) this.vertList[6] = this.interpolate(6); // 6 7
        if ((edgeInfo & 128) > 0) this.vertList[7] = this.interpolate(7); // 4 7
        if ((edgeInfo & 256) > 0) this.vertList[8] = this.interpolate(8); // 0 4
        if ((edgeInfo & 512) > 0) this.vertList[9] = this.interpolate(9); // 1 5
        if ((edgeInfo & 1024) > 0) this.vertList[10] = this.interpolate(10); // 2 6
        if ((edgeInfo & 2048) > 0) this.vertList[11] = this.interpolate(11); // 3 7

        let triInfo = TriTable[tableIndex];
        for (let t = 0; t < triInfo.length; t += 3) {
            this.triangleCount++;
            ChunkedArray.add3(this.triangleBuffer, this.vertList[triInfo[t]], this.vertList[triInfo[t + 1]], this.vertList[triInfo[t + 2]]);
        }
    }
}