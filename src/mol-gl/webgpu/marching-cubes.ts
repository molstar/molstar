/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { fillSerial } from '../../mol-util/array';
import { MarchingCubesParams } from '../../mol-geo/util/marching-cubes/algorithm';
import { CubeEdges, CubeVertices, TriTable } from '../../mol-geo/util/marching-cubes/tables';
import { RuntimeContext } from '../../mol-task';
import { WebGPUContext } from './context';
import { GPUBufferUsage, GPUMapMode } from './compat';

export interface WebGPUPackedScalarField {
    readonly buffer: GPUBuffer
    readonly count: number
    readonly dimensions: readonly number[]
    dispose?(): void
}

const corners = CubeVertices.map(p => `vec3i(${p.i}, ${p.j}, ${p.k})`).join(',');
const edges = CubeEdges.map(e => `vec2u(${CubeVertices.indexOf(e.a)}, ${CubeVertices.indexOf(e.b)})`).join(',');
const triangles = TriTable.flatMap(row => Array.from({ length: 16 }, (_, i) => row[i] ?? -1)).join(',');
const CellFloats = 184;
export const marchingCubesShader = /* wgsl */ `
struct Parameters { dimensions: vec4u, strides: vec4u, extent: vec4u, batch: vec4u, shape: vec4f, periodicity: vec4u, region: vec4u };
struct Vertex { position: vec4f, normal: vec4f, group: vec4f };
struct Cell { info: vec4u, vertices: array<Vertex, 15> };
@group(0) @binding(0) var<uniform> p: Parameters;
@group(0) @binding(1) var<storage, read> field: array<f32>;
@group(0) @binding(2) var<storage, read> ids: array<f32>;
@group(0) @binding(3) var<storage, read_write> cells: array<Cell>;
const corners = array<vec3i, 8>(${corners});
const edges = array<vec2u, 12>(${edges});
const triangles = array<i32, 4096>(${triangles});
fn offset(c: vec3i) -> u32 { let coordinate = vec3u(c) % p.periodicity.xyz; return coordinate.x * p.strides.x + coordinate.y * p.strides.y + coordinate.z * p.strides.z; }
fn sample(c: vec3i) -> f32 { return field[offset(clamp(c, vec3i(0), vec3i(p.dimensions.xyz) - vec3i(1)))]; }
fn gradient(c: vec3i) -> vec3f {
    return vec3f(sample(c - vec3i(1,0,0)) - sample(c + vec3i(1,0,0)), sample(c - vec3i(0,1,0)) - sample(c + vec3i(0,1,0)), sample(c - vec3i(0,0,1)) - sample(c + vec3i(0,0,1)));
}
@compute @workgroup_size(64)
fn compute(@builtin(global_invocation_id) invocation: vec3u) {
    let local = invocation.x + invocation.y * p.batch.z * 64u;
    if (local >= p.batch.y) { return; }
    let index = p.batch.x + local;
    let coordinate = vec3i(vec3u(index % p.extent.x, (index / p.extent.x) % p.extent.y, index / (p.extent.x * p.extent.y))) + vec3i(p.region.xyz);
    var table = 0u;
    for (var i = 0u; i < 8u; i++) { let value = sample(coordinate + corners[i]); if (value < p.shape.x || (value == p.shape.x && p.shape.y > 0.0)) { table |= 1u << i; } }
    var count = 0u;
    for (var triangle = 0u; triangle < 15u; triangle += 3u) {
        if (triangles[table * 16u + triangle] < 0) { break; }
        var vertices: array<Vertex, 3>;
        var ignored = false;
        for (var v = 0u; v < 3u; v++) {
            let order = select(v, 2u - v, p.shape.x < 0.0);
            let edge = edges[u32(triangles[table * 16u + triangle + order])];
            let a = coordinate + corners[edge.x]; let b = coordinate + corners[edge.y];
            let va = sample(a); let vb = sample(b);
            let t = (p.shape.x - va) / (va - vb);
            var group = 0.0;
            if (p.strides.w == 1u) { group = f32(offset(coordinate)); }
            if (p.strides.w == 2u) {
                let ga = ids[offset(a)]; let gb = ids[offset(b)];
                group = select(gb, ga, t < 0.5);
                if (group == -1.0) { group = select(ga, gb, t < 0.5); }
                if (group == -2.0) { ignored = true; }
            }
            let normal = gradient(a) + t * (gradient(a) - gradient(b));
            vertices[v] = Vertex(vec4f(vec3f(a) + t * vec3f(a - b), 1.0), vec4f(normal * select(-1.0, 1.0, p.shape.x >= 0.0), 0.0), vec4f(group, 0.0, 0.0, 0.0));
        }
        if (!ignored) {
            for (var v = 0u; v < 3u; v++) { cells[local].vertices[count + v] = vertices[v]; }
            count += 3u;
        }
    }
    cells[local].info = vec4u(count, 0u, 0u, 0u);
}
`;

// Gaussian density compute stores scalar values and closest-atom IDs as one
// vec2f per voxel. This variant lets marching cubes consume that buffer
// directly without a GPU-to-CPU round trip or a second upload.
const packedMarchingCubesShader = marchingCubesShader
    .replace('field: array<f32>', 'field: array<vec2f>')
    .replace('ids: array<f32>', 'ids: array<vec2f>')
    .replace('@group(0) @binding(2) var<storage, read> ids: array<vec2f>;\n', '')
    .replace('field[offset(clamp(c, vec3i(0), vec3i(p.dimensions.xyz) - vec3i(1)))]', 'field[offset(clamp(c, vec3i(0), vec3i(p.dimensions.xyz) - vec3i(1)))].x')
    .replace(/ids\[offset\((a|b)\)\]/g, 'field[offset($1)].y');

/** Stable prefix scan and compaction; no atomics or unordered vertex appends. */
export const marchingCubesCompactionShader = /* wgsl */ `
struct Parameters { dimensions: vec4u, strides: vec4u, extent: vec4u, batch: vec4u, shape: vec4f, periodicity: vec4u, region: vec4u };
struct Vertex { position: vec4f, normal: vec4f, group: vec4f };
struct Cell { info: vec4u, vertices: array<Vertex, 15> };
@group(0) @binding(0) var<uniform> p: Parameters;
@group(0) @binding(1) var<storage, read> cells: array<Cell>;
@group(0) @binding(2) var<storage, read_write> offsets: array<u32>;
@group(0) @binding(3) var<storage, read_write> blocks: array<u32>;
@group(0) @binding(4) var<storage, read_write> compacted: array<Vertex>;
var<workgroup> prefix: array<u32, 64>;
@compute @workgroup_size(64)
fn scan(@builtin(workgroup_id) workgroup: vec3u, @builtin(local_invocation_index) lane: u32) {
    let block = workgroup.x + workgroup.y * p.batch.z;
    if (block >= (p.batch.y + 63u) / 64u) { return; }
    let index = block * 64u + lane;
    var count = 0u;
    if (index < p.batch.y) { count = cells[index].info.x; }
    prefix[lane] = count;
    workgroupBarrier();
    for (var step = 1u; step < 64u; step *= 2u) {
        var previous = 0u;
        if (lane >= step) { previous = prefix[lane - step]; }
        workgroupBarrier();
        prefix[lane] += previous;
        workgroupBarrier();
    }
    if (index < p.batch.y) { offsets[index] = prefix[lane] - count; }
    if (lane == 63u) { blocks[block] = prefix[lane]; }
}
@compute @workgroup_size(1)
fn scanBlocks() {
    var total = 0u;
    for (var i = 0u; i < (p.batch.y + 63u) / 64u; i++) {
        let count = blocks[i]; blocks[i] = total; total += count;
    }
    offsets[p.batch.y] = total;
}
@compute @workgroup_size(64)
fn compact(@builtin(global_invocation_id) invocation: vec3u) {
    let index = invocation.x + invocation.y * p.batch.z * 64u;
    if (index >= p.batch.y) { return; }
    let start = offsets[index] + blocks[index / 64u];
    for (var i = 0u; i < cells[index].info.x; i++) { compacted[start + i] = cells[index].vertices[i]; }
}
`;

interface CompactionPipelines { scan: GPUComputePipeline, scanBlocks: GPUComputePipeline, compact: GPUComputePipeline, layout: GPUBindGroupLayout }
const compactionPipelines = new WeakMap<GPUDevice, Promise<CompactionPipelines>>();
function getCompactionPipelines(device: GPUDevice) {
    let result = compactionPipelines.get(device);
    if (!result) {
        result = (async () => {
            const module = device.createShaderModule({ label: 'molstar-marching-cubes-compaction', code: marchingCubesCompactionShader });
            const errors = (await module.getCompilationInfo()).messages.filter(m => m.type === 'error');
            if (errors.length) throw new Error(`WebGPU marching cubes compaction: ${errors.map(m => m.message).join('; ')}`);
            const layout = device.createBindGroupLayout({ entries: [0, 1, 2, 3, 4].map(binding => ({ binding, visibility: 4, buffer: { type: binding === 0 ? 'uniform' : binding === 1 ? 'read-only-storage' : 'storage' } })) });
            const pipelineLayout = device.createPipelineLayout({ bindGroupLayouts: [layout] });
            const [scan, scanBlocks, compact] = await Promise.all(['scan', 'scanBlocks', 'compact'].map(entryPoint => device.createComputePipelineAsync({ layout: pipelineLayout, compute: { module, entryPoint } })));
            return { scan, scanBlocks, compact, layout };
        })();
        compactionPipelines.set(device, result);
    }
    return result;
}

const pipelines = new WeakMap<GPUDevice, { normal?: Promise<GPUComputePipeline>, packed?: Promise<GPUComputePipeline> }>();
async function pipeline(device: GPUDevice, packed = false) {
    const entries = pipelines.get(device) ?? {};
    const key = packed ? 'packed' : 'normal';
    let result = entries[key];
    if (!result) {
        result = (async () => {
            const module = device.createShaderModule({ label: packed ? 'molstar-marching-cubes-packed' : 'molstar-marching-cubes', code: packed ? packedMarchingCubesShader : marchingCubesShader });
            const errors = (await module.getCompilationInfo()).messages.filter(m => m.type === 'error');
            if (errors.length) throw new Error(`WebGPU marching cubes shader: ${errors.map(e => e.message).join('; ')}`);
            return device.createComputePipelineAsync({ layout: 'auto', compute: { module, entryPoint: 'compute' } });
        })();
        entries[key] = result; pipelines.set(device, entries);
    }
    return result;
}

/** GPU-generated triangle soup in grid coordinates, compacted before bounded readback.
 * With retainGPU, only counts are downloaded; the caller owns the returned GPU buffer.
 */
export async function computeMarchingCubesWebGPU(ctx: RuntimeContext, context: WebGPUContext, params: MarchingCubesParams, options: { maxCellsPerBatch?: number, groupMode?: 'zero' | 'cell' | 'id', periodicDimensions?: ReadonlyArray<number>, retainGPU?: boolean } = {}) {
    const packed = (params as MarchingCubesParams & { webgpuField?: WebGPUPackedScalarField }).webgpuField;
    const dimensions = params.scalarField.space.dimensions;
    if (dimensions.length !== 3 || !dimensions.every(d => Number.isSafeInteger(d) && d >= 2) || !Number.isFinite(Math.fround(params.isoLevel))) throw new Error('WebGPU marching cubes requires a finite iso-level and a 3D grid with at least two samples per axis.');
    const bottomLeft = params.bottomLeft ?? [0, 0, 0], topRight = params.topRight ?? dimensions;
    if (bottomLeft.length !== 3 || topRight.length !== 3 || !bottomLeft.every((d, i) => Number.isSafeInteger(d) && d >= 0 && Number.isSafeInteger(topRight[i]) && d <= topRight[i] && topRight[i] <= dimensions[i])) throw new Error('WebGPU marching cubes region must lie within the scalar grid.');
    const periodicity = options.periodicDimensions ?? dimensions;
    if (periodicity.length !== 3 || !periodicity.every((d, i) => Number.isSafeInteger(d) && d > 0 && (options.periodicDimensions ? dimensions[i] === d + 1 : dimensions[i] === d))) throw new Error('WebGPU marching cubes periodic dimensions must match the wrapped grid.');
    if (params.scalarField.data.length !== periodicity[0] * periodicity[1] * periodicity[2]) throw new Error('WebGPU marching cubes scalar data must match its grid dimensions.');
    // Tensor views (including floodfill) can change samples without changing
    // their backing array. Read through the accessor in original storage order.
    if (packed && (packed.count !== params.scalarField.data.length || packed.dimensions.some((d, i) => d !== dimensions[i]))) throw new Error('WebGPU packed scalar field does not match the marching-cubes grid.');
    const field = packed ? undefined : new Float32Array(params.scalarField.data.length), coordinate = [0, 0, 0];
    if (field) {
        for (let i = 0; i < field.length; i++) {
            params.scalarField.space.getCoords(i, coordinate);
            field[i] = params.scalarField.space.get(params.scalarField.data, ...coordinate);
        }
        if (!field.every(Number.isFinite)) throw new Error('WebGPU marching cubes requires finite scalar samples.');
    }
    const strides = [0, 1, 2].map(a => { const coordinate = [0, 0, 0]; coordinate[a] = 1; return params.scalarField.space.dataOffset(...coordinate); });
    const mode = packed ? 'id' : options.groupMode ?? (params.idField ? 'id' : 'zero');
    const ids = new Float32Array(mode === 'id' && !packed ? params.scalarField.data.length : 1);
    if (mode === 'id' && !packed) {
        if (!params.idField || !params.idField.space.dimensions.every((d, i) => d === dimensions[i])) throw new Error('WebGPU marching cubes IDs must match the scalar grid.');
        const coordinate = [0, 0, 0];
        for (let i = 0; i < ids.length; i++) { params.scalarField.space.getCoords(i, coordinate); ids[i] = params.idField.space.get(params.idField.data, ...coordinate); }
        if (!ids.every(Number.isFinite)) throw new Error('WebGPU marching cubes requires finite IDs.');
    }
    const extent = dimensions.map((_, i) => Math.max(0, topRight[i] - bottomLeft[i] - 1)), count = extent[0] * extent[1] * extent[2];
    const maxBytes = Math.min(context.device.limits.maxStorageBufferBindingSize, context.device.limits.maxBufferSize);
    if (!Number.isSafeInteger(count) || count > 0xffffffff || (!packed && field!.byteLength > maxBytes) || (!packed && ids.byteLength > maxBytes) || (packed && packed.buffer.size > maxBytes)) throw new Error('WebGPU marching cubes grid exceeds the device storage limit.');
    const requested = options.maxCellsPerBatch ?? 4096;
    if (!Number.isSafeInteger(requested) || requested < 1) throw new Error('WebGPU marching cubes batch size must be a positive integer.');
    if (count === 0) return { positions: new Float32Array(0), normals: new Float32Array(0), groups: new Float32Array(0), vertexCount: 0, readbackBytes: 0, gpuBuffer: undefined as GPUBuffer | undefined };
    const batchSize = Math.min(count, requested, Math.floor(maxBytes / (CellFloats * 4)));
    const { device } = context, computePipeline = await pipeline(device, !!packed), compaction = await getCompactionPipelines(device), buffers: GPUBuffer[] = [];
    const chunks: { positions: Float32Array, normals: Float32Array, groups: Float32Array }[] = [];
    const gpuChunks: GPUBuffer[] = [];
    let vertexCount = 0, readbackBytes = 0;
    let gpuBuffer: GPUBuffer | undefined;
    try {
        const scalars = packed ? packed.buffer : context.createBuffer(field!, GPUBufferUsage.STORAGE); if (!packed) buffers.push(scalars);
        const groups = packed ? packed.buffer : context.createBuffer(ids, GPUBufferUsage.STORAGE); if (!packed) buffers.push(groups);
        const uniform = device.createBuffer({ size: 112, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST }); buffers.push(uniform);
        const output = device.createBuffer({ size: batchSize * CellFloats * 4, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC }); buffers.push(output);
        const offsets = device.createBuffer({ size: (batchSize + 1) * 4, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC }); buffers.push(offsets);
        const blocks = device.createBuffer({ size: Math.ceil(batchSize / 64) * 4, usage: GPUBufferUsage.STORAGE }); buffers.push(blocks);
        const compacted = device.createBuffer({ size: batchSize * 15 * 48, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC }); buffers.push(compacted);
        const readback = options.retainGPU ? undefined : device.createBuffer({ size: compacted.size, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ }); if (readback) buffers.push(readback);
        const countReadback = device.createBuffer({ size: 4, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ }); buffers.push(countReadback);
        const compactionGroup = device.createBindGroup({ layout: compaction.layout, entries: [uniform, output, offsets, blocks, compacted].map((buffer, binding) => ({ binding, resource: { buffer } })) });
        const bindGroupEntries = packed
            ? [{ binding: 0, resource: { buffer: uniform } }, { binding: 1, resource: { buffer: scalars } }, { binding: 3, resource: { buffer: output } }]
            : [uniform, scalars, groups, output].map((buffer, binding) => ({ binding, resource: { buffer } }));
        const bindGroup = device.createBindGroup({ layout: computePipeline.getBindGroupLayout(0), entries: bindGroupEntries });
        const data = new Float32Array(28), uints = new Uint32Array(data.buffer);
        uints.set(dimensions); uints.set([...strides, mode === 'id' ? 2 : mode === 'cell' ? 1 : 0], 4); uints.set(extent, 8); data[16] = params.isoLevel; data[17] = params.isoLevel - Math.fround(params.isoLevel); uints.set(periodicity, 20); uints.set(bottomLeft, 24);
        for (let base = 0; base < count; base += batchSize) {
            const size = Math.min(batchSize, count - base), x = Math.min(Math.ceil(size / 64), device.limits.maxComputeWorkgroupsPerDimension);
            uints.set([base, size, x, 0], 12); device.queue.writeBuffer(uniform, 0, data);
            const encoder = device.createCommandEncoder(); const pass = encoder.beginComputePass();
            pass.setPipeline(computePipeline); pass.setBindGroup(0, bindGroup); pass.dispatchWorkgroups(x, Math.ceil(size / (x * 64))); pass.end();
            for (const stage of ['scan', 'scanBlocks', 'compact'] as const) {
                const compactPass = encoder.beginComputePass(); compactPass.setPipeline(compaction[stage]); compactPass.setBindGroup(0, compactionGroup);
                if (stage === 'scanBlocks') compactPass.dispatchWorkgroups(1);
                else compactPass.dispatchWorkgroups(x, Math.ceil(size / (x * 64)));
                compactPass.end();
            }
            encoder.copyBufferToBuffer(offsets, size * 4, countReadback, 0, 4);
            device.queue.submit([encoder.finish()]); context.stats.computeDispatches += 4; context.stats.marchingCubesDispatches += 4;
            await countReadback.mapAsync(GPUMapMode.READ);
            const total = new Uint32Array(countReadback.getMappedRange())[0]; countReadback.unmap(); readbackBytes += 4;
            if (total > size * 15 || total % 3 !== 0) throw new Error('WebGPU marching cubes produced an invalid compacted vertex count.');
            const positions = new Float32Array(options.retainGPU ? 0 : total * 3), normals = new Float32Array(options.retainGPU ? 0 : total * 3), groups = new Float32Array(options.retainGPU ? 0 : total);
            if (total && options.retainGPU) {
                if ((vertexCount + total) * 80 > maxBytes) throw new Error('Native surface geometry exceeds the device storage buffer limit.');
                const chunk = device.createBuffer({ size: total * 48, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC });
                gpuChunks.push(chunk);
                const copy = device.createCommandEncoder(); copy.copyBufferToBuffer(compacted, 0, chunk, 0, total * 48); device.queue.submit([copy.finish()]);
            } else if (total && readback) {
                const copy = device.createCommandEncoder(); copy.copyBufferToBuffer(compacted, 0, readback, 0, total * 48); device.queue.submit([copy.finish()]);
                await readback.mapAsync(GPUMapMode.READ, 0, total * 48);
                const values = new Float32Array(readback.getMappedRange(0, total * 48));
                for (let i = 0; i < total; i++) {
                    const offset = i * 12;
                    positions.set(values.subarray(offset, offset + 3), i * 3); normals.set(values.subarray(offset + 4, offset + 7), i * 3); groups[i] = values[offset + 8];
                }
                readback.unmap(); readbackBytes += total * 48;
            }
            chunks.push({ positions, normals, groups }); vertexCount += total;
            if (ctx.shouldUpdate) await ctx.update({ message: 'Extracting surface on WebGPU', current: base + size, max: count });
        }
        if (options.retainGPU) {
            gpuBuffer = device.createBuffer({ label: 'molstar-native-surface-cells', size: Math.max(48, vertexCount * 48), usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
            const encoder = device.createCommandEncoder(); let offset = 0;
            for (const chunk of gpuChunks) { encoder.copyBufferToBuffer(chunk, 0, gpuBuffer, offset, chunk.size); offset += chunk.size; }
            device.queue.submit([encoder.finish()]);
        }
    } catch (error) { gpuBuffer?.destroy(); throw error; } finally { for (const buffer of [...buffers, ...gpuChunks]) buffer.destroy(); }
    if (options.retainGPU) return { positions: new Float32Array(0), normals: new Float32Array(0), groups: new Float32Array(0), vertexCount, readbackBytes, gpuBuffer };
    const positions = new Float32Array(vertexCount * 3), normals = new Float32Array(vertexCount * 3), groups = new Float32Array(vertexCount);
    let offset = 0; for (const chunk of chunks) { positions.set(chunk.positions, offset * 3); normals.set(chunk.normals, offset * 3); groups.set(chunk.groups, offset); offset += chunk.groups.length; }
    return { positions, normals, groups, vertexCount, readbackBytes, gpuBuffer };
}

/** Shared representation bridge while direct GPU geometry consumption is implemented. */
export async function computeMarchingCubesMeshWebGPU(ctx: RuntimeContext, context: WebGPUContext, params: MarchingCubesParams, mesh?: Mesh, groupMode?: 'zero' | 'cell' | 'id', periodicDimensions?: ReadonlyArray<number>) {
    const result = await computeMarchingCubesWebGPU(ctx, context, params, { groupMode, periodicDimensions });
    return Mesh.create(result.positions, fillSerial(new Uint32Array(result.vertexCount)), result.normals, result.groups, result.vertexCount, result.vertexCount / 3, mesh);
}
