/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from '@molstar/core/task';
import { Tensor } from '@molstar/core/math/linear-algebra';
import type { Mesh } from '../../geometry/mesh/mesh.js';
import { type Index, EdgeIdInfo, CubeEdges, EdgeTable, TriTable } from './tables.js';
import { defaults } from '@molstar/core/util';
import { anchorMap } from '@molstar/core/util/map-anchor';
import { type MarchingCubesBuilder, MarchingCubesMeshBuilder, MarchingCubesLinesBuilder } from './builder.js';
import type { Lines } from '../../geometry/lines/lines.js';

/**
 * The parameters required by the algorithm.
 */
export interface MarchingCubesParams {
  isoLevel: number;
  scalarField: Tensor;
  bottomLeft?: ReadonlyArray<number>;
  topRight?: ReadonlyArray<number>;
  idField?: Tensor;
}

interface MarchingCubesInputParams extends MarchingCubesParams {
  bottomLeft: ReadonlyArray<number>;
  topRight: ReadonlyArray<number>;
}

function getInputParams(params: MarchingCubesParams): MarchingCubesInputParams {
  return {
    ...params,
    bottomLeft: defaults(params.bottomLeft, [0, 0, 0] as ReadonlyArray<number>),
    topRight: defaults(params.topRight, params.scalarField.space.dimensions),
  };
}

function getExtent(inputParams: MarchingCubesInputParams) {
  return {
    dX: inputParams.topRight[0] - inputParams.bottomLeft[0],
    dY: inputParams.topRight[1] - inputParams.bottomLeft[1],
    dZ: inputParams.topRight[2] - inputParams.bottomLeft[2],
  };
}

export function computeMarchingCubesMesh(params: MarchingCubesParams, mesh?: Mesh) {
  return Task.create('Marching Cubes Mesh', async (ctx) => {
    const inputParams = getInputParams(params);
    const { dX, dY, dZ } = getExtent(inputParams);
    // TODO should it be configurable? Scalar fields can produce meshes with vastly different densities.
    const vertexChunkSize = Math.min(262144, Math.max((dX * dY * dZ) / 32, 1024));
    const builder = MarchingCubesMeshBuilder(vertexChunkSize, mesh);
    await new MarchingCubesComputation(ctx, builder, inputParams).run();
    return builder.get();
  });
}

export function computeMarchingCubesLines(params: MarchingCubesParams, lines?: Lines) {
  return Task.create('Marching Cubes Lines', async (ctx) => {
    const inputParams = getInputParams(params);
    const { dX, dY, dZ } = getExtent(inputParams);
    // TODO should it be configurable? Scalar fields can produce meshes with vastly different densities.
    const vertexChunkSize = Math.min(262144, Math.max((dX * dY * dZ) / 32, 1024));
    const builder = MarchingCubesLinesBuilder(vertexChunkSize, lines);
    await new MarchingCubesComputation(ctx, builder, inputParams).run();
    return builder.get();
  });
}

class MarchingCubesComputation {
  private size: number;
  private sliceSize: number;
  // set here, not first in `doSlices`, so that a computation has its final map when made (`anchorMap` below)
  private edgeFilter = 0;

  private minX = 0;
  private minY = 0;
  private minZ = 0;
  private maxX = 0;
  private maxY = 0;
  private maxZ = 0;
  private state: MarchingCubesState;

  private async doSlices() {
    let done = 0;

    this.edgeFilter = 15;
    for (let k = this.minZ; k < this.maxZ; k++, this.edgeFilter &= ~4) {
      this.slice(k);

      done += this.sliceSize;
      if (this.ctx.shouldUpdate) {
        await this.ctx.update({ message: 'Computing surface...', current: done, max: this.size });
      }
    }
  }

  private slice(k: number) {
    this.edgeFilter |= 2;
    for (let j = this.minY; j < this.maxY; j++, this.edgeFilter &= ~2) {
      // bit 1 is set for the first cell of the row only (`processRow` clears it)
      this.state.processRow(this.minX, this.maxX, j, k, this.edgeFilter | 1);
    }
    this.state.clearEdgeVertexIndexSlice(k);
  }

  async run() {
    await this.doSlices();
  }

  constructor(
    private ctx: RuntimeContext,
    builder: MarchingCubesBuilder<any>,
    params: MarchingCubesInputParams,
  ) {
    this.state = new MarchingCubesState(builder, params);
    this.minX = params.bottomLeft[0];
    this.minY = params.bottomLeft[1];
    this.minZ = params.bottomLeft[2];
    this.maxX = params.topRight[0] - 1;
    this.maxY = params.topRight[1] - 1;
    this.maxZ = params.topRight[2] - 1;

    this.size = (this.maxX - this.minX) * (this.maxY - this.minY) * (this.maxZ - this.minZ);
    this.sliceSize = (this.maxX - this.minX) * (this.maxY - this.minY);
  }
}

class MarchingCubesState {
  nX: number;
  nY: number;
  nZ: number;
  isoLevel: number;
  scalarFieldGet: Tensor.Space['get'];
  scalarField: Tensor.Data;
  idFieldGet?: Tensor.Space['get'];
  idField?: Tensor.Data;

  /**
   * Whether the scalar field is read directly at strided offsets (`Tensor.strides`) instead of through `get`.
   * Both read the same element, so the output is the same; the direct read avoids a closure call per value, and the
   * call site going polymorphic when fields of several layouts (or floodfilled ones) are meshed in one session.
   */
  strided: boolean;
  sx = 0;
  sy = 0;
  sz = 0;
  /** Offsets of the cube corners 1 to 7 (`CubeVertices`) from corner 0, when strided */
  c1 = 0;
  c2 = 0;
  c3 = 0;
  c4 = 0;
  c5 = 0;
  c6 = 0;
  c7 = 0;
  /** Offsets of the two ends of each cube edge (`CubeEdges`) from corner 0, when strided */
  edgeA: number[] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
  edgeB: number[] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

  idStrided: boolean;
  idSx = 0;
  idSy = 0;
  idSz = 0;

  // two layers of vertex indices. Each vertex has 3 edges associated.
  verticesOnEdges: Int32Array;
  vertList: number[] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
  i: number = 0;
  j: number = 0;
  k: number = 0;
  /** The data offset of the current cell, when strided */
  o: number = 0;

  private get3dOffsetFromEdgeInfo(index: Index) {
    return this.nX * (((this.k + index.k) % 2) * this.nY + this.j + index.j) + this.i + index.i;
  }

  /**
   * This clears the "vertex index buffer" for the slice that will not be accessed anymore.
   */
  clearEdgeVertexIndexSlice(k: number) {
    // clear either the top or bottom half of the buffer...
    const start = k % 2 === 0 ? 0 : 3 * this.nX * this.nY;
    const end = k % 2 === 0 ? 3 * this.nX * this.nY : this.verticesOnEdges.length;
    this.verticesOnEdges.fill(0, start, end);
  }

  private interpolate(edgeNum: number) {
    const info = EdgeIdInfo[edgeNum];
    const edgeId = 3 * this.get3dOffsetFromEdgeInfo(info) + info.e;

    const ret = this.verticesOnEdges[edgeId];
    if (ret > 0) return ret - 1;

    const sf = this.scalarField;
    const sfg = this.scalarFieldGet;
    const strided = this.strided;

    const edge = CubeEdges[edgeNum];
    const a = edge.a,
      b = edge.b;
    const li = a.i + this.i,
      lj = a.j + this.j,
      lk = a.k + this.k;
    const hi = b.i + this.i,
      hj = b.j + this.j,
      hk = b.k + this.k;
    // data offsets of the edge ends, when strided
    const lo = this.o + this.edgeA[edgeNum];
    const ho = this.o + this.edgeB[edgeNum];
    const v0 = strided ? sf[lo] : sfg(sf, li, lj, lk);
    const v1 = strided ? sf[ho] : sfg(sf, hi, hj, hk);
    const t = (this.isoLevel - v0) / (v0 - v1);

    if (this.idField) {
      const idf = this.idField;
      let u: number, v: number;
      if (this.idStrided) {
        const { idSx, idSy, idSz } = this;
        u = idf[li * idSx + lj * idSy + lk * idSz];
        v = idf[hi * idSx + hj * idSy + hk * idSz];
      } else {
        u = this.idFieldGet!(idf, li, lj, lk);
        v = this.idFieldGet!(idf, hi, hj, hk);
      }
      let a = t < 0.5 ? u : v;
      // -1 means 'no id', check if the other cell has an id
      if (a === -1) a = t < 0.5 ? v : u;
      // -2 means 'ignore this cell'
      if (a === -2) return -1;
      this.builder.addGroup(a);
    } else {
      this.builder.addGroup(0);
    }

    const id = this.builder.addVertex(li + t * (li - hi), lj + t * (lj - hj), lk + t * (lk - hk));
    this.verticesOnEdges[edgeId] = id + 1;

    // TODO cache scalarField differences for slices
    // TODO make calculation optional
    let n0x: number, n0y: number, n0z: number, n1x: number, n1y: number, n1z: number;
    if (strided) {
      // central differences, clamped at the faces as `Math.max(0, i - 1)` and `Math.min(nX - 1, i + 1)` below
      const { sx, sy, sz } = this;
      const mX = this.nX - 1,
        mY = this.nY - 1,
        mZ = this.nZ - 1;
      n0x = sf[li > 0 ? lo - sx : lo] - sf[li < mX ? lo + sx : lo];
      n0y = sf[lj > 0 ? lo - sy : lo] - sf[lj < mY ? lo + sy : lo];
      n0z = sf[lk > 0 ? lo - sz : lo] - sf[lk < mZ ? lo + sz : lo];
      n1x = sf[hi > 0 ? ho - sx : ho] - sf[hi < mX ? ho + sx : ho];
      n1y = sf[hj > 0 ? ho - sy : ho] - sf[hj < mY ? ho + sy : ho];
      n1z = sf[hk > 0 ? ho - sz : ho] - sf[hk < mZ ? ho + sz : ho];
    } else {
      n0x = sfg(sf, Math.max(0, li - 1), lj, lk) - sfg(sf, Math.min(this.nX - 1, li + 1), lj, lk);
      n0y = sfg(sf, li, Math.max(0, lj - 1), lk) - sfg(sf, li, Math.min(this.nY - 1, lj + 1), lk);
      n0z = sfg(sf, li, lj, Math.max(0, lk - 1)) - sfg(sf, li, lj, Math.min(this.nZ - 1, lk + 1));

      n1x = sfg(sf, Math.max(0, hi - 1), hj, hk) - sfg(sf, Math.min(this.nX - 1, hi + 1), hj, hk);
      n1y = sfg(sf, hi, Math.max(0, hj - 1), hk) - sfg(sf, hi, Math.min(this.nY - 1, hj + 1), hk);
      n1z = sfg(sf, hi, hj, Math.max(0, hk - 1)) - sfg(sf, hi, hj, Math.min(this.nZ - 1, hk + 1));
    }

    const nx = n0x + t * (n0x - n1x);
    const ny = n0y + t * (n0y - n1y);
    const nz = n0z + t * (n0z - n1z);

    // ensure normal-direction is the same for negative and positive iso-levels
    if (this.isoLevel >= 0) {
      this.builder.addNormal(nx, ny, nz);
    } else {
      this.builder.addNormal(-nx, -ny, -nz);
    }

    return id;
  }

  constructor(
    private builder: MarchingCubesBuilder<any>,
    params: MarchingCubesInputParams,
  ) {
    const dims = params.scalarField.space.dimensions;
    this.nX = dims[0];
    this.nY = dims[1];
    this.nZ = dims[2];
    this.isoLevel = params.isoLevel;
    this.scalarFieldGet = params.scalarField.space.get;
    this.scalarField = params.scalarField.data;

    const strides = Tensor.strides(params.scalarField.space);
    // the clamped central differences in `interpolate` match `get`'s only for cells inside the grid
    const { bottomLeft: bl, topRight: tr } = params;
    const inside = bl[0] >= 0 && bl[1] >= 0 && bl[2] >= 0 && tr[0] <= dims[0] && tr[1] <= dims[1] && tr[2] <= dims[2];
    this.strided = !!strides && strides.length === 3 && inside;
    if (this.strided) {
      const [sx, sy, sz] = strides!;
      this.sx = sx;
      this.sy = sy;
      this.sz = sz;
      this.c1 = sx;
      this.c2 = sx + sy;
      this.c3 = sy;
      this.c4 = sz;
      this.c5 = sx + sz;
      this.c6 = sx + sy + sz;
      this.c7 = sy + sz;
      for (let e = 0; e < 12; e++) {
        const { a, b } = CubeEdges[e];
        this.edgeA[e] = a.i * sx + a.j * sy + a.k * sz;
        this.edgeB[e] = b.i * sx + b.j * sy + b.k * sz;
      }
    }

    this.idField = undefined;
    this.idFieldGet = undefined;
    this.idStrided = false;
    if (params.idField) {
      this.idField = params.idField.data;
      this.idFieldGet = params.idField.space.get;
      const idStrides = Tensor.strides(params.idField.space);
      if (idStrides && idStrides.length === 3) {
        this.idStrided = true;
        [this.idSx, this.idSy, this.idSz] = idStrides;
      }
    }

    // two layers of vertex indices. Each vertex has 3 edges associated.
    this.verticesOnEdges = new Int32Array(3 * this.nX * this.nY * 2);
  }

  /**
   * The cells `(minX..maxX - 1, j, k)`. Each cell reads only the corners of its far face (i + 1) and takes those of its
   * near face from the previous cell; when strided, the data offset is carried along i.
   */
  processRow(minX: number, maxX: number, j: number, k: number, edgeFilter: number) {
    if (minX >= maxX) return;
    if (this.strided) this.stridedRow(minX, maxX, j, k, edgeFilter);
    else this.getRow(minX, maxX, j, k, edgeFilter);
  }

  private stridedRow(minX: number, maxX: number, j: number, k: number, edgeFilter: number) {
    const iso = this.isoLevel;
    const sf = this.scalarField;
    const { sx, c1, c2, c3, c4, c5, c6, c7 } = this;
    let o = minX * sx + j * this.sy + k * this.sz;

    // corners 0, 3, 4, 7 (`CubeVertices`) of the first cell
    let near = 0;
    if (sf[o] < iso) near |= 1;
    if (sf[o + c3] < iso) near |= 8;
    if (sf[o + c4] < iso) near |= 16;
    if (sf[o + c7] < iso) near |= 128;

    for (let i = minX; i < maxX; i++, o += sx, edgeFilter &= ~1) {
      // corners 1, 2, 5, 6
      let far = 0;
      if (sf[o + c1] < iso) far |= 2;
      if (sf[o + c2] < iso) far |= 4;
      if (sf[o + c5] < iso) far |= 32;
      if (sf[o + c6] < iso) far |= 64;

      const tableIndex = near | far;
      // the far face is the next cell's near face: corners 1 → 0, 2 → 3, 5 → 4, 6 → 7
      near = ((far & 2) >> 1) | ((far & 4) << 1) | ((far & 32) >> 1) | ((far & 64) << 1);

      if (tableIndex === 0 || tableIndex === 255) continue;
      this.processCell(i, j, k, o, tableIndex, edgeFilter);
    }
  }

  /** As `stridedRow`, through `get` */
  private getRow(minX: number, maxX: number, j: number, k: number, edgeFilter: number) {
    const iso = this.isoLevel;
    const sf = this.scalarField;
    const get = this.scalarFieldGet;

    let near = 0;
    if (get(sf, minX, j, k) < iso) near |= 1;
    if (get(sf, minX, j + 1, k) < iso) near |= 8;
    if (get(sf, minX, j, k + 1) < iso) near |= 16;
    if (get(sf, minX, j + 1, k + 1) < iso) near |= 128;

    for (let i = minX; i < maxX; i++, edgeFilter &= ~1) {
      const i1 = i + 1;
      let far = 0;
      if (get(sf, i1, j, k) < iso) far |= 2;
      if (get(sf, i1, j + 1, k) < iso) far |= 4;
      if (get(sf, i1, j, k + 1) < iso) far |= 32;
      if (get(sf, i1, j + 1, k + 1) < iso) far |= 64;

      const tableIndex = near | far;
      near = ((far & 2) >> 1) | ((far & 4) << 1) | ((far & 32) >> 1) | ((far & 64) << 1);

      if (tableIndex === 0 || tableIndex === 255) continue;
      this.processCell(i, j, k, 0, tableIndex, edgeFilter);
    }
  }

  private processCell(i: number, j: number, k: number, o: number, tableIndex: number, edgeFilter: number) {
    this.i = i;
    this.j = j;
    this.k = k;
    this.o = o;
    const edgeInfo = EdgeTable[tableIndex];
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

    const triInfo = TriTable[tableIndex];
    for (let t = 0; t < triInfo.length; t += 3) {
      const l = triInfo[t],
        m = triInfo[t + 1],
        n = triInfo[t + 2];
      // ensure winding-order is the same for negative and positive iso-levels
      if (this.isoLevel >= 0) {
        this.builder.addTriangle(this.vertList, l, m, n, edgeFilter);
      } else {
        this.builder.addTriangle(this.vertList, n, m, l, edgeFilter);
      }
    }
  }
}

// keeps the maps of a computation, its state and the mesh builder, which the optimized code of `slice`, `processCell`,
// `interpolate` & co. embeds (`anchorMap`); made, not run, so that it adds no feedback (perf-audit-3 08/P15)
{
  const space = Tensor.Space([2, 2, 2], [0, 1, 2], Float32Array);
  const params = getInputParams({ isoLevel: 0.5, scalarField: Tensor.create(space, space.create()) });
  anchorMap(new MarchingCubesComputation(RuntimeContext.Synchronous, MarchingCubesMeshBuilder(1), params));
}
