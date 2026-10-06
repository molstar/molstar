/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { OrderedSet } from '../../mol-data/int/ordered-set';
import { RuntimeContext } from '../../mol-task';
import { Box3D } from '../../mol-math/geometry/primitives/box3d';
import type { PositionData } from '../../mol-math/geometry/common';
import type { GaussianDensityData, GaussianDensityProps } from '../../mol-math/geometry/gaussian-density';
import { Mat4, Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { WebGPUContext } from './context';
import { GPUBufferUsage, GPUMapMode } from './compat';

export interface WebGPUGaussianDensityBuffer {
    readonly buffer: GPUBuffer
    readonly count: number
    readonly dimensions: readonly number[]
    dispose(): void
}

export const gaussianDensityShader = /* wgsl */ `
struct Parameters {
    dimensions: vec4u,
    originScale: vec4f,
    batch: vec4u,
    shape: vec4f,
    bins: vec4u,
};
struct Atom { positionRadius: vec4f, info: vec4f };
@group(0) @binding(0) var<uniform> parameters: Parameters;
@group(0) @binding(1) var<storage, read> atoms: array<Atom>;
@group(0) @binding(2) var<storage, read_write> density: array<vec2f>;
@group(0) @binding(3) var<storage, read> binRanges: array<vec2u>;
@group(0) @binding(4) var<storage, read> binAtoms: array<u32>;
@compute @workgroup_size(64)
fn compute(@builtin(global_invocation_id) invocation: vec3u) {
    let local = invocation.x + invocation.y * parameters.batch.z * 64u;
    if (local >= parameters.batch.y) { return; }
    let index = parameters.batch.x + local;
    let dimensions = parameters.dimensions.xyz;
    let coordinate = vec3u(index / (dimensions.y * dimensions.z), (index / dimensions.z) % dimensions.y, index % dimensions.z);
    let point = parameters.originScale.xyz + vec3f(coordinate) * parameters.originScale.w;
    var total = 0.0; var strongest = 0.0; var id = -1.0;
    let cell = coordinate / parameters.batch.w;
    let range = binRanges[(cell.x * parameters.bins.y + cell.y) * parameters.bins.z + cell.z];
    for (var i = 0u; i < range.y; i++) {
        let atom = atoms[binAtoms[range.x + i]];
        let delta = point - atom.positionRadius.xyz;
        let distanceSquared = dot(delta, delta);
        let radiusSquared = atom.positionRadius.w * atom.positionRadius.w;
        if (radiusSquared == 0.0) { continue; }
        if (distanceSquared <= 4.0 * radiusSquared) {
            let contribution = exp(-parameters.shape.x * distanceSquared / radiusSquared);
            total += contribution;
            // Strict comparison keeps the first atom in equal-density ties.
            if (contribution > strongest) { strongest = contribution; id = atom.info.x; }
        }
    }
    density[local] = vec2f(total, id);
}
`;

const pipelines = new WeakMap<GPUDevice, Promise<GPUComputePipeline>>();
function getPipeline(device: GPUDevice) {
    let pipeline = pipelines.get(device);
    if (!pipeline) {
        pipeline = (async () => {
            const module = device.createShaderModule({ label: 'molstar-gaussian-density', code: gaussianDensityShader });
            const errors = (await module.getCompilationInfo()).messages.filter(m => m.type === 'error');
            if (errors.length) throw new Error(`WebGPU Gaussian density shader: ${errors.map(e => e.message).join('; ')}`);
            return device.createComputePipelineAsync({ label: 'molstar-gaussian-density', layout: 'auto', compute: { module, entryPoint: 'compute' } });
        })();
        pipelines.set(device, pipeline);
    }
    return pipeline;
}

/** Conservative voxel bins. Lists retain input order, preserving sums and tie-breaking. */
export function createGaussianDensityBins(atoms: Float32Array, count: number, dimensions: ArrayLike<number>, origin: ArrayLike<number>, resolution: number, maxBytes: number, enabled = true) {
    if (!Number.isFinite(Math.fround(resolution)) || Math.fround(resolution) <= 0) throw new Error('WebGPU Gaussian resolution must be representable as a positive f32 value.');
    const limit = Math.min(0xffffffff, Math.floor(maxBytes / 4));
    if (limit < Math.max(2, count)) throw new Error('WebGPU Gaussian bins exceed the device storage buffer limit.');
    const largestDimension = Math.max(dimensions[0], dimensions[1], dimensions[2]);
    let maxRadius = 0;
    for (let i = 0; i < count; i++) maxRadius = Math.max(maxRadius, atoms[i * 8 + 3]);
    let span = enabled ? Math.min(largestDimension, Math.max(1, Math.ceil(2 * maxRadius / resolution))) : largestDimension;
    for (;;) {
        const dims = [0, 1, 2].map(a => Math.ceil(dimensions[a] / span));
        const cells = dims[0] * dims[1] * dims[2];
        if (cells * 2 > limit) { span = Math.min(largestDimension, span * 2); continue; }
        const ranges = new Uint32Array(cells * 2);
        // Include rounding in both world-to-grid conversion and the shader's f32
        // multiply/add/subtract. Padding extra voxels is cheap and conservative.
        const bounds = new Float64Array(count * 6);
        for (let i = 0; i < count; i++) {
            for (let a = 0; a < 3; a++) {
                const start = Math.fround(origin[a]), step = Math.fround(resolution), center = atoms[i * 8 + a];
                const extent = Math.max(Math.abs(start), Math.abs(start + dimensions[a] * step), Math.abs(center), 1);
                const pad = 2 + extent * 2 ** -20 / step;
                const coordinate = (center - start) / step, cutoff = 2 * atoms[i * 8 + 3] / step;
                bounds[i * 6 + a] = !enabled ? 0 : Math.max(0, Math.min(dims[a], Math.floor((coordinate - cutoff - pad) / span)));
                bounds[i * 6 + a + 3] = !enabled ? dims[a] - 1 : Math.max(-1, Math.min(dims[a] - 1, Math.floor((coordinate + cutoff + pad) / span)));
            }
        }
        const visit = (i: number, fn: (cell: number) => void) => {
            if (atoms[i * 8 + 3] === 0) return;
            const b = i * 6;
            for (let x = bounds[b]; x <= bounds[b + 3]; x++) {
                for (let y = bounds[b + 1]; y <= bounds[b + 4]; y++) {
                    for (let z = bounds[b + 2]; z <= bounds[b + 5]; z++) fn((x * dims[1] + y) * dims[2] + z);
                }
            }
        };
        let entries = 0;
        for (let i = 0; i < count && entries <= limit; i++) visit(i, cell => { ranges[cell * 2 + 1]++; entries++; });
        if (entries > limit) { span = Math.min(largestDimension, span * 2); continue; }
        const indices = new Uint32Array(Math.max(1, entries)), cursors = new Uint32Array(cells);
        let offset = 0, evaluations = 0;
        for (let x = 0; x < dims[0]; x++) for (let y = 0; y < dims[1]; y++) for (let z = 0; z < dims[2]; z++) {
            const cell = (x * dims[1] + y) * dims[2] + z;
            ranges[cell * 2] = offset; cursors[cell] = offset;
            const length = ranges[cell * 2 + 1]; offset += length;
            evaluations += length * Math.min(span, dimensions[0] - x * span) * Math.min(span, dimensions[1] - y * span) * Math.min(span, dimensions[2] - z * span);
        }
        for (let i = 0; i < count; i++) visit(i, cell => { indices[cursors[cell]++] = i; });
        return { ranges, indices, dimensions: dims, span, evaluations };
    }
}

/** Density and closest-atom grids, in the same z-fastest tensor order as the CPU path. */
export async function GaussianDensityWebGPU(ctx: RuntimeContext, context: WebGPUContext, position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, options: { maxVoxelsPerBatch?: number, spatialBins?: boolean } = {}): Promise<GaussianDensityData> {
    const { resolution, radiusOffset, smoothness } = props;
    if (!Number.isFinite(resolution) || resolution <= 0 || !Number.isFinite(smoothness) || smoothness <= 0) throw new Error('WebGPU Gaussian density requires positive finite resolution and smoothness.');
    const n = OrderedSet.size(position.indices);
    const atoms = new Float32Array(Math.max(1, n) * 8);
    let maxRadius = 0;
    for (let i = 0; i < n; i++) {
        const index = OrderedSet.getAt(position.indices, i), r = radius(index) + radiusOffset;
        if (!Number.isFinite(r) || r < 0) throw new Error('WebGPU Gaussian density requires nonnegative finite atom radii.');
        if (![position.x[index], position.y[index], position.z[index], position.id?.[i] ?? i].every(Number.isFinite)) throw new Error('WebGPU Gaussian density requires finite atom coordinates and IDs.');
        maxRadius = Math.max(maxRadius, r);
        atoms.set([position.x[index], position.y[index], position.z[index], r, position.id?.[i] ?? i], i * 8);
    }
    if (!atoms.every(Number.isFinite)) throw new Error('WebGPU Gaussian atom data must be representable as finite f32 values.');
    const pad = maxRadius * 2 + resolution;
    const expandedBox = Box3D.expand(Box3D(), box, Vec3.create(pad, pad, pad));
    if (!expandedBox.min.every(v => Number.isFinite(Math.fround(v)))) throw new Error('WebGPU Gaussian grid origin must be representable as finite f32 values.');
    const scaledBox = Box3D.scale(Box3D(), expandedBox, 1 / resolution);
    const dimensions = Vec3.ceil(Vec3(), Box3D.size(Vec3(), scaledBox));
    const count = dimensions[0] * dimensions[1] * dimensions[2];
    if (!dimensions.every(d => Number.isSafeInteger(d) && d > 0) || !Number.isSafeInteger(count) || count > 0xffffffff) throw new Error('WebGPU Gaussian density grid exceeds the addressable voxel range.');
    const { device } = context;
    if (atoms.byteLength > device.limits.maxStorageBufferBindingSize || atoms.byteLength > device.limits.maxBufferSize) throw new Error('WebGPU Gaussian atom data exceeds the device storage buffer limit.');
    const bins = createGaussianDensityBins(atoms, n, dimensions, expandedBox.min, resolution, Math.min(device.limits.maxStorageBufferBindingSize, device.limits.maxBufferSize), options.spatialBins !== false);
    const requestedBatch = options.maxVoxelsPerBatch ?? 262144;
    if (!Number.isSafeInteger(requestedBatch) || requestedBatch < 1) throw new Error('WebGPU Gaussian density batch size must be a positive integer.');
    const batchSize = Math.min(count, requestedBatch, Math.floor(device.limits.maxStorageBufferBindingSize / 8), Math.floor(device.limits.maxBufferSize / 8));
    const space = Tensor.Space(dimensions, [0, 1, 2], Float32Array);
    const data = space.create(), ids = space.create();
    const pipeline = await getPipeline(device);
    const buffers: GPUBuffer[] = [];
    let retainedDensity: GPUBuffer | undefined;
    try {
        const atomBuffer = context.createBuffer(atoms, GPUBufferUsage.STORAGE, 'molstar-gaussian-atoms'); buffers.push(atomBuffer);
        const uniform = device.createBuffer({ size: 80, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST }); buffers.push(uniform);
        const output = device.createBuffer({ size: batchSize * 8, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC }); buffers.push(output);
        const readback = device.createBuffer({ size: batchSize * 8, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST }); buffers.push(readback);
        const retainDensity = count * 8 <= Math.min(device.limits.maxStorageBufferBindingSize, device.limits.maxBufferSize);
        const densityGrid = retainDensity ? device.createBuffer({ label: 'molstar-gaussian-density-grid', size: Math.max(8, count * 8), usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST }) : undefined;
        if (densityGrid) buffers.push(densityGrid);
        const ranges = context.createBuffer(bins.ranges, GPUBufferUsage.STORAGE, 'molstar-gaussian-bin-ranges'); buffers.push(ranges);
        const indices = context.createBuffer(bins.indices, GPUBufferUsage.STORAGE, 'molstar-gaussian-bin-atoms'); buffers.push(indices);
        const bindGroup = device.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [
            { binding: 0, resource: { buffer: uniform } }, { binding: 1, resource: { buffer: atomBuffer } }, { binding: 2, resource: { buffer: output } },
            { binding: 3, resource: { buffer: ranges } }, { binding: 4, resource: { buffer: indices } },
        ] });
        const parameters = new Float32Array(20), uints = new Uint32Array(parameters.buffer);
        uints.set([...dimensions, n]); parameters.set([...expandedBox.min, resolution], 4); parameters[12] = smoothness; uints.set(bins.dimensions, 16);
        for (let base = 0; base < count; base += batchSize) {
            const size = Math.min(batchSize, count - base);
            const workgroupsX = Math.min(Math.ceil(size / 64), device.limits.maxComputeWorkgroupsPerDimension);
            const workgroupsY = Math.ceil(size / (workgroupsX * 64));
            uints.set([base, size, workgroupsX, bins.span], 8);
            device.queue.writeBuffer(uniform, 0, parameters);
            const encoder = device.createCommandEncoder({ label: 'molstar-gaussian-density' });
            const pass = encoder.beginComputePass(); pass.setPipeline(pipeline); pass.setBindGroup(0, bindGroup);
            pass.dispatchWorkgroups(workgroupsX, workgroupsY); pass.end();
            encoder.copyBufferToBuffer(output, 0, readback, 0, size * 8);
            if (densityGrid) encoder.copyBufferToBuffer(output, 0, densityGrid, base * 8, size * 8);
            device.queue.submit([encoder.finish()]);
            context.stats.computeDispatches++;
            await readback.mapAsync(GPUMapMode.READ, 0, size * 8);
            const values = new Float32Array(readback.getMappedRange(0, size * 8));
            for (let j = 0; j < size; j++) { data[base + j] = values[j * 2]; ids[base + j] = values[j * 2 + 1]; }
            readback.unmap();
            if (ctx.shouldUpdate) await ctx.update({ message: 'computing density grid on WebGPU', current: base + size, max: count });
        }
        retainedDensity = densityGrid;
    } finally { for (const buffer of buffers) if (buffer !== retainedDensity) buffer.destroy(); }
    const transform = Mat4.fromScaling(Mat4(), Vec3.create(resolution, resolution, resolution));
    Mat4.setTranslation(transform, expandedBox.min);
    const result = { field: Tensor.create(space, data), idField: Tensor.create(space, ids), transform, radiusFactor: 1, resolution, maxRadius } as GaussianDensityData & { webgpuDensity?: WebGPUGaussianDensityBuffer };
    if (retainedDensity) {
        result.webgpuDensity = { buffer: retainedDensity, count, dimensions: [...dimensions], dispose: () => retainedDensity?.destroy() };
    }
    return result;
}
