/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { TextureMesh } from '../../mol-geo/geometry/texture-mesh/texture-mesh';
import { MarchingCubesParams } from '../../mol-geo/util/marching-cubes/algorithm';
import { Sphere3D } from '../../mol-math/geometry';
import { Mat4 } from '../../mol-math/linear-algebra';
import { RuntimeContext } from '../../mol-task';
import { packIntToRGBArray } from '../../mol-util/number-packing';
import { GPUBufferUsage, GPUMapMode, GPUTextureUsage } from './compat';
import { WebGPUContext } from './context';
import type { WebGPUGeometry } from './geometry';
import { computeMarchingCubesWebGPU } from './marching-cubes';
import { WebGPUTextureData } from './texture-data';

const shader = /* wgsl */ `
struct Parameters { transform: mat4x4f, normalTransform: mat4x4f, shape: vec4u };
struct Vertex { position: vec4f, normal: vec4f, group: vec4f };
struct GeometryVertex { position: vec4f, normal: vec4f, end: vec4f, mapping: vec4f, info: vec4f };
@group(0) @binding(0) var<uniform> p: Parameters;
@group(0) @binding(1) var<storage, read> input: array<Vertex>;
@group(0) @binding(2) var<storage, read_write> output: array<GeometryVertex>;
@group(0) @binding(3) var positions: texture_storage_2d<rgba32float, write>;
@group(0) @binding(4) var normals: texture_storage_2d<rgba32float, write>;
@group(0) @binding(5) var groups: texture_storage_2d<rgba8unorm, write>;
@compute @workgroup_size(64)
fn compute(@builtin(global_invocation_id) invocation: vec3u) {
    let index = invocation.x + invocation.y * p.shape.w * 64u;
    if (index >= p.shape.y * p.shape.z) { return; }
    let pixel = vec2i(i32(index % p.shape.y), i32(index / p.shape.y));
    var position = vec4f(0.0); var normal = vec4f(0.0); var packedGroup = vec4f(0.0);
    if (index < p.shape.x) {
        position = p.transform * input[index].position;
        var n = (p.normalTransform * vec4f(input[index].normal.xyz, 0.0)).xyz;
        let lengthSquared = dot(n, n);
        if (lengthSquared > 0.0) { n *= inverseSqrt(lengthSquared); }
        normal = vec4f(n, 0.0);
        let group = input[index].group.x;
        let encoded = u32(max(group, 0.0)) + 1u;
        packedGroup = vec4f(f32((encoded >> 16u) & 255u), f32((encoded >> 8u) & 255u), f32(encoded & 255u), 0.0) / 255.0;
        output[index] = GeometryVertex(position, normal, vec4f(0.0), vec4f(0.0), vec4f(group, f32(index), 0.0, 0.0));
    }
    textureStore(positions, pixel, position); textureStore(normals, pixel, normal); textureStore(groups, pixel, packedGroup);
}
`;

const pipelines = new WeakMap<GPUDevice, Promise<GPUComputePipeline>>();
async function getPipeline(device: GPUDevice) {
    let pipeline = pipelines.get(device);
    if (!pipeline) {
        pipeline = (async () => {
            const module = device.createShaderModule({ label: 'molstar-native-texture-mesh', code: shader });
            const errors = (await module.getCompilationInfo()).messages.filter(m => m.type === 'error');
            if (errors.length) throw new Error(`Native texture mesh shader: ${errors.map(m => m.message).join('; ')}`);
            return device.createComputePipelineAsync({ layout: 'auto', compute: { module, entryPoint: 'compute' } });
        })();
        pipelines.set(device, pipeline);
    }
    return pipeline;
}

/** GPU-generated storage vertices and textures; the CPU mirror serves synchronous theme/location APIs. */
export class WebGPUTextureMeshGeometry {
    private disposed = false;
    private readonly versions: number[];
    readonly geometry: WebGPUGeometry;
    constructor(readonly device: GPUDevice, readonly buffer: GPUBuffer, vertices: Float32Array, readonly position: WebGPUTextureData, readonly normal: WebGPUTextureData, readonly group: WebGPUTextureData) {
        const vertexCount = vertices.length / 20;
        this.geometry = { vertices, indices: Uint32Array.from({ length: vertexCount }, (_, i) => i), vertexCount, kind: 0, native: { device, buffer } };
        this.versions = [position, normal, group].map(t => t.version);
    }
    matches(position: unknown, normal: unknown, group: unknown, vertexCount: number) {
        return !this.disposed && this.position === position && this.normal === normal && this.group === group && vertexCount === this.geometry.vertexCount && [this.position, this.normal, this.group].every((t, i) => t.version === this.versions[i]);
    }
    destroy() {
        if (this.disposed) return;
        this.disposed = true; this.buffer.destroy();
        this.position.destroy(); this.normal.destroy(); this.group.destroy();
    }
}

/** Extract, compact and transform on the GPU. Rendering reuses the compute output without a vertex upload. */
export async function computeMarchingCubesTextureMeshWebGPU(ctx: RuntimeContext, context: WebGPUContext, params: MarchingCubesParams, transform: Mat4, groupCount: number, boundingSphere: Sphere3D, previous?: TextureMesh, options: { maxCellsPerBatch?: number } = {}) {
    if (!Number.isSafeInteger(groupCount) || groupCount < 0 || groupCount > 16777215) throw new Error('Native texture mesh group count exceeds RGB packing.');
    if (!transform.every(Number.isFinite)) throw new Error('Native texture mesh transform must be finite.');
    const inverse = Mat4();
    if (!Mat4.tryInvert(inverse, transform)) throw new Error('Native texture mesh transform must be invertible.');
    const normalTransform = Mat4.transpose(Mat4(), inverse);
    const { device } = context;
    const result = await computeMarchingCubesWebGPU(ctx, context, params, { ...options, retainGPU: true });
    const input = result.gpuBuffer ?? context.createBuffer(new Float32Array(12), GPUBufferUsage.STORAGE);
    const count = result.vertexCount, dimensionLimit = device.limits.maxTextureDimension2D;
    const width = Math.max(1, Math.min(dimensionLimit, Math.ceil(Math.sqrt(count)))), height = Math.max(1, Math.ceil(count / width));
    const buffers: GPUBuffer[] = [input], textures: GPUTexture[] = [];
    let geometry: WebGPUTextureMeshGeometry | undefined;
    try {
        if (height > dimensionLimit || count * 80 > device.limits.maxStorageBufferBindingSize) throw new Error('Native texture mesh exceeds device storage/texture limits.');
        const pipeline = await getPipeline(device);
        const vertexBuffer = device.createBuffer({ label: 'molstar-native-surface-vertices', size: Math.max(80, count * 80), usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC }); buffers.push(vertexBuffer);
        const readback = device.createBuffer({ size: vertexBuffer.size, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST }); buffers.push(readback);
        const data = new Float32Array(36), workgroups = Math.min(Math.ceil(width * height / 64), device.limits.maxComputeWorkgroupsPerDimension);
        data.set(transform); data.set(normalTransform, 16); new Uint32Array(data.buffer).set([count, width, height, workgroups], 32);
        const uniform = context.createBuffer(data, GPUBufferUsage.UNIFORM); buffers.push(uniform);
        for (const format of ['rgba32float', 'rgba32float', 'rgba8unorm'] as const) textures.push(device.createTexture({ label: 'molstar-native-surface-texture', size: [width, height], format, usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC }));
        const bindGroup = device.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [
            { binding: 0, resource: { buffer: uniform } }, { binding: 1, resource: { buffer: input } }, { binding: 2, resource: { buffer: vertexBuffer } },
            ...textures.map((texture, i) => ({ binding: i + 3, resource: texture.createView() })),
        ] });
        const encoder = device.createCommandEncoder({ label: 'molstar-native-surface-pack' });
        const pass = encoder.beginComputePass(); pass.setPipeline(pipeline); pass.setBindGroup(0, bindGroup); pass.dispatchWorkgroups(workgroups, Math.ceil(width * height / (workgroups * 64))); pass.end();
        encoder.copyBufferToBuffer(vertexBuffer, 0, readback, 0, vertexBuffer.size); device.queue.submit([encoder.finish()]); context.stats.computeDispatches++;
        await readback.mapAsync(GPUMapMode.READ);
        const vertices = new Float32Array(readback.getMappedRange()).slice(0, count * 20); readback.unmap();
        const positions = new Float32Array(width * height * 4), normals = new Float32Array(positions.length), groups = new Uint8Array(positions.length);
        for (let i = 0; i < count; i++) {
            const offset = i * 20;
            if (!vertices.subarray(offset, offset + 8).every(Number.isFinite)) throw new Error('Native texture mesh geometry must be finite.');
            const group = vertices[offset + 16];
            if (!Number.isSafeInteger(group) || group < 0 || group >= groupCount) throw new Error('Native texture mesh group lies outside its location domain.');
            positions.set(vertices.subarray(offset, offset + 4), i * 4); normals.set(vertices.subarray(offset + 4, offset + 8), i * 4); packIntToRGBArray(group, groups, i * 4);
        }
        const native = [new WebGPUTextureData(), new WebGPUTextureData(), new WebGPUTextureData()];
        for (let i = 0; i < 3; i++) native[i].loadGPU(device, textures[i], { width, height, depth: 1, array: [positions, normals, groups][i] });
        geometry = new WebGPUTextureMeshGeometry(device, vertexBuffer, vertices, native[0], native[1], native[2]);
        const old = previous?.meta.webgpuGeometry;
        const surface = TextureMesh.create(count, groupCount, native[0], native[2], native[1], boundingSphere, previous);
        surface.meta.webgpuGeometry = geometry;
        surface.meta.webgpuReadbackBytes = result.readbackBytes + vertexBuffer.size;
        if (old instanceof WebGPUTextureMeshGeometry) old.destroy();
        return surface;
    } catch (error) {
        geometry?.destroy(); for (const texture of textures) texture.destroy(); throw error;
    } finally {
        for (const buffer of buffers) if (buffer !== geometry?.buffer) buffer.destroy();
    }
}
