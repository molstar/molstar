/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { ColorSmoothingInput, ColorSmoothingOptions } from '../../mol-geo/geometry/mesh/color-smoothing';
import { TextureMeshValues } from '../renderable/texture-mesh';
import { unpackRGBToInt } from '../../mol-util/number-packing';
import { MeshValues } from '../renderable/mesh';
import { Box3D } from '../../mol-math/geometry';
import { Vec2, Vec3, Vec4 } from '../../mol-math/linear-algebra';
import { getVolumeTexture2dLayout } from '../../mol-repr/volume/util';
import { ValueCell } from '../../mol-util';
import { GPUBufferUsage, GPUMapMode, GPUTextureUsage } from './compat';
import { WebGPUContext } from './context';
import { WebGPUTextureData } from './texture-data';
import { createGaussianDensityBins } from './gaussian-density';

const smoothingShader = /* wgsl */ `
struct Parameters { dimensions: vec4u, batch: vec4u, bins: vec4u, texture: vec4u };
struct Sample { position: vec4f, color: vec4f };
@group(0) @binding(0) var<uniform> params: Parameters;
@group(0) @binding(1) var<storage, read> samples: array<Sample>;
@group(0) @binding(2) var<storage, read> ranges: array<vec2u>;
@group(0) @binding(3) var<storage, read> indices: array<u32>;
fn gridColor(coordinate: vec3u) -> vec4u {
    let cell = coordinate / params.batch.w;
    let range = ranges[(cell.x * params.bins.y + cell.y) * params.bins.z + cell.z];
    var sum = vec4f(0.0); var weight = 0.0;
    for (var i = 0u; i < range.y; i++) {
        let sample = samples[indices[range.x + i]];
        let d = distance(vec3f(coordinate), sample.position.xyz);
        if (d > 2.0) { continue; }
        let w = 2.0 - d;
        sum += sample.color * w; weight += w;
    }
    var color = vec4u(0u);
    if (weight > 0.0) { color = vec4u(clamp(floor(sum / weight + 0.5), vec4f(0.0), vec4f(255.0))); }
    if (params.dimensions.w == 3u) { color.a = 255u; }
    return color;
}
`;
const shader = smoothingShader + /* wgsl */ `
@group(0) @binding(4) var<storage, read_write> output: array<u32>;
@compute @workgroup_size(64)
fn compute(@builtin(global_invocation_id) invocation: vec3u) {
    let local = invocation.x + invocation.y * params.batch.z * 64u;
    if (local >= params.batch.y) { return; }
    let index = params.batch.x + local;
    let dim = params.dimensions.xyz;
    let coordinate = vec3u(index / (dim.y * dim.z), (index / dim.z) % dim.y, index % dim.z);
    let color = gridColor(coordinate);
    output[local] = color.r | (color.g << 8u) | (color.b << 16u) | (color.a << 24u);
}
`;
const textureShader = smoothingShader + /* wgsl */ `
@group(0) @binding(4) var output: texture_storage_2d<rgba8unorm, write>;
@compute @workgroup_size(64)
fn compute(@builtin(global_invocation_id) invocation: vec3u) {
    let index = invocation.x + invocation.y * params.batch.z * 64u;
    if (index >= params.batch.y) { return; }
    let pixel = vec2u(index % params.texture.x, index / params.texture.x);
    let dim = params.dimensions.xyz;
    let columns = params.texture.x / dim.x;
    let coordinate = vec3u(pixel.x % dim.x, pixel.y % dim.y, (pixel.y / dim.y) * columns + pixel.x / dim.x);
    var color = vec4u(0u);
    if (params.dimensions.w == 3u) { color.a = 255u; }
    if (coordinate.z < dim.z) { color = gridColor(coordinate); }
    textureStore(output, vec2i(pixel), vec4f(color) / 255.0);
}
`;
const pipelines = new WeakMap<GPUDevice, Promise<GPUComputePipeline>>();
function getPipeline(device: GPUDevice) {
    let pipeline = pipelines.get(device);
    if (!pipeline) {
        pipeline = (async () => {
            const module = device.createShaderModule({ label: 'molstar-color-smoothing', code: shader });
            const errors = (await module.getCompilationInfo()).messages.filter(m => m.type === 'error');
            if (errors.length) throw new Error(`WebGPU smoothing shader: ${errors.map(e => e.message).join('; ')}`);
            return device.createComputePipelineAsync({ layout: 'auto', compute: { module, entryPoint: 'compute' } });
        })();
        pipelines.set(device, pipeline);
    }
    return pipeline;
}

function prepareSmoothing(context: WebGPUContext, input: ColorSmoothingInput, options: ColorSmoothingOptions) {
    const { resolution, stride } = options;
    if (!Number.isFinite(resolution) || resolution <= 0 || !Number.isSafeInteger(stride) || stride < 1) throw new Error('WebGPU smoothing requires positive resolution and stride.');
    const instanceType = input.colorType === 'groupInstance';
    const instanceCount = instanceType ? input.instanceCount : 1;
    if (![input.vertexCount, input.instanceCount, input.groupCount].every(n => Number.isSafeInteger(n) && n >= 0) || input.positionBuffer.length < input.vertexCount * 3 || input.groupBuffer.length < input.vertexCount || (instanceType && (input.transformBuffer.length < instanceCount * 16 || input.instanceBuffer.length < instanceCount))) throw new Error('Invalid WebGPU smoothing geometry.');
    const box = Box3D.fromSphere3D(Box3D(), instanceType ? input.boundingSphere : input.invariantBoundingSphere);
    const pad = 1 + resolution;
    const expanded = Box3D.expand(Box3D(), box, Vec3.create(pad, pad, pad));
    const scale = 1 / resolution;
    const gridDim = Vec3.ceil(Vec3(), Box3D.size(Vec3(), Box3D.scale(Box3D(), expanded, scale)));
    Vec3.add(gridDim, gridDim, Vec3.create(2, 2, 2));
    const count = gridDim[0] * gridDim[1] * gridDim[2];
    if (!gridDim.every(d => Number.isSafeInteger(d) && d > 0) || !Number.isSafeInteger(count) || count > 0xffffffff || ![...expanded.min, scale].every(v => Number.isFinite(Math.fround(v)))) throw new Error('Invalid WebGPU smoothing grid.');
    const { width, height } = getVolumeTexture2dLayout(gridDim);
    const sampleCount = instanceCount * Math.ceil(input.vertexCount / stride);
    const { device } = context;
    if (width > device.limits.maxTextureDimension2D || height > device.limits.maxTextureDimension2D) throw new Error('WebGPU smoothing grid exceeds the device texture dimension limit.');
    const maxBytes = Math.min(device.limits.maxStorageBufferBindingSize, device.limits.maxBufferSize);
    if (!Number.isSafeInteger(sampleCount) || Math.max(1, sampleCount) * 32 > maxBytes) throw new Error('WebGPU smoothing samples exceed the device buffer limit.');
    const samples = new Float32Array(Math.max(1, sampleCount) * 8), point = Vec3();
    let offset = 0;
    for (let i = 0; i < instanceCount; i++) {
        const instance = instanceType ? input.instanceBuffer[i] : 0;
        if (instanceType && (!Number.isInteger(instance) || instance < 0 || instance >= input.instanceCount)) throw new Error('Invalid WebGPU smoothing instance ID.');
        for (let j = 0; j < input.vertexCount; j += stride) {
            const group = input.groupBuffer[j], color = (instance * input.groupCount + group) * input.itemSize;
            if (!Number.isInteger(group) || group < 0 || group >= input.groupCount || color + input.itemSize > input.colorData.array.length) throw new Error('Invalid WebGPU smoothing group or theme data.');
            Vec3.fromArray(point, input.positionBuffer, j * 3);
            if (instanceType) Vec3.transformMat4Offset(point, point, input.transformBuffer, 0, 0, i * 16);
            Vec3.sub(point, point, expanded.min); Vec3.scale(point, point, scale);
            samples.set([...point, 1], offset);
            const colors = input.colorData.array;
            if (input.itemSize === 1) samples[offset + 7] = colors[color];
            else for (let k = 0; k < input.itemSize; k++) samples[offset + 4 + k] = colors[color + k];
            offset += 8;
        }
    }
    if (!samples.every(Number.isFinite)) throw new Error('WebGPU smoothing samples must be finite f32 values.');
    // The bin builder uses radius * 2. Unit radii give the smoothing kernel's
    // two-grid-unit support, including conservative floating-point padding.
    const bins = createGaussianDensityBins(samples, sampleCount, gridDim, [0, 0, 0], 1, maxBytes);
    return { device, maxBytes, width, height, count, gridDim, samples, bins, gridTransform: Vec4.create(expanded.min[0], expanded.min[1], expanded.min[2], scale), type: instanceType ? 'volumeInstance' as const : 'volume' as const };
}

/** Native weighted grid accumulation. Ordered spatial bins preserve sample order. */
export async function calcMeshColorSmoothingWebGPU(context: WebGPUContext, input: ColorSmoothingInput, options: ColorSmoothingOptions, texture = new WebGPUTextureData(), batchLimit = 262144) {
    if (!Number.isSafeInteger(batchLimit) || batchLimit < 1) throw new Error('WebGPU smoothing batch size must be positive.');
    const { device, maxBytes, width, height, count, gridDim, samples, bins, gridTransform, type } = prepareSmoothing(context, input, options);
    const batchSize = Math.min(count, batchLimit, Math.floor(maxBytes / 4), device.limits.maxComputeWorkgroupsPerDimension ** 2 * 64);
    const grid = new Uint8Array(width * height * 4);
    if (input.itemSize === 3) for (let i = 3; i < grid.length; i += 4) grid[i] = 255;
    const pipeline = await getPipeline(device), buffers: GPUBuffer[] = [];
    try {
        const inputs = [samples, bins.ranges, bins.indices].map(a => { const b = context.createBuffer(a, GPUBufferUsage.STORAGE); buffers.push(b); return b; });
        const uniform = device.createBuffer({ size: 64, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST }); buffers.push(uniform);
        const output = device.createBuffer({ size: batchSize * 4, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC }); buffers.push(output);
        const readback = device.createBuffer({ size: batchSize * 4, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST }); buffers.push(readback);
        const bindGroup = device.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [uniform, ...inputs, output].map((buffer, binding) => ({ binding, resource: { buffer } })) });
        const params = new Uint32Array(16); params.set([...gridDim, input.itemSize]); params.set(bins.dimensions, 8);
        for (let base = 0; base < count; base += batchSize) {
            const size = Math.min(batchSize, count - base), x = Math.min(Math.ceil(size / 64), device.limits.maxComputeWorkgroupsPerDimension);
            params.set([base, size, x, bins.span], 4); device.queue.writeBuffer(uniform, 0, params);
            const encoder = device.createCommandEncoder({ label: 'molstar-color-smoothing' });
            const pass = encoder.beginComputePass(); pass.setPipeline(pipeline); pass.setBindGroup(0, bindGroup); pass.dispatchWorkgroups(x, Math.ceil(size / (x * 64))); pass.end();
            encoder.copyBufferToBuffer(output, 0, readback, 0, size * 4); device.queue.submit([encoder.finish()]); context.stats.computeDispatches++;
            await readback.mapAsync(GPUMapMode.READ, 0, size * 4);
            const bytes = new Uint8Array(readback.getMappedRange(0, size * 4));
            for (let i = 0; i < size; i++) {
                const index = base + i, z = index % gridDim[2], y = Math.floor(index / gridDim[2]) % gridDim[1], gx = Math.floor(index / (gridDim[1] * gridDim[2]));
                const target = ((Math.floor(z * gridDim[0] / width) * gridDim[1] + y) * width + (z * gridDim[0]) % width + gx) * 4;
                grid.set(bytes.subarray(i * 4, i * 4 + 4), target);
            }
            readback.unmap();
        }
    } finally { buffers.forEach(b => b.destroy()); }
    texture.load({ array: grid, width, height });
    return { kind: 'volume' as const, texture, gridTexDim: Vec2.create(width, height), gridDim, gridTransform, type };
}

function smoothingGeometry(values: MeshValues | TextureMeshValues) {
    if ('aPosition' in values) return { positionBuffer: values.aPosition.ref.value, groupBuffer: values.aGroup.ref.value };
    const position = values.tPosition.ref.value, group = values.tGroup.ref.value;
    if (!(position instanceof WebGPUTextureData) || !(group instanceof WebGPUTextureData)) throw new Error('Native texture-mesh smoothing requires WebGPU texture data.');
    const positions = position.data.array, groups = group.data.array, vertexCount = values.uVertexCount.ref.value;
    if (!(positions instanceof Float32Array) || !(groups instanceof Float32Array || groups instanceof Uint8Array)) throw new Error('Unsupported native smoothing geometry.');
    const positionBuffer = new Float32Array(vertexCount * 3), groupBuffer = new Float32Array(vertexCount), scale = groups instanceof Uint8Array ? 1 : 255;
    for (let i = 0; i < vertexCount; i++) {
        positionBuffer.set(positions.subarray(i * 4, i * 4 + 3), i * 3);
        groupBuffer[i] = unpackRGBToInt(Math.round(groups[i * 4] * scale), Math.round(groups[i * 4 + 1] * scale), Math.round(groups[i * 4 + 2] * scale));
    }
    return { positionBuffer, groupBuffer };
}

export async function applyMeshColorSmoothingWebGPU(context: WebGPUContext, values: MeshValues | TextureMeshValues, options: ColorSmoothingOptions, texture?: WebGPUTextureData) {
    const colorType = values.dColorType.ref.value;
    if (colorType !== 'group' && colorType !== 'groupInstance') return;
    const data = calcMeshColorSmoothingTextureWebGPU(context, {
        vertexCount: values.uVertexCount.ref.value, instanceCount: values.uInstanceCount.ref.value, groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value, instanceBuffer: values.aInstance.ref.value, ...smoothingGeometry(values),
        colorData: values.tColor.ref.value, colorType, boundingSphere: values.boundingSphere.ref.value, invariantBoundingSphere: values.invariantBoundingSphere.ref.value, itemSize: 3,
    }, options, texture);
    ValueCell.updateIfChanged(values.dColorType, data.type); ValueCell.update(values.tColorGrid, data.texture);
    ValueCell.update(values.uColorTexDim, data.gridTexDim); ValueCell.update(values.uColorGridDim, data.gridDim); ValueCell.update(values.uColorGridTransform, data.gridTransform);
}

const texturePipelines = new WeakMap<GPUDevice, GPUComputePipeline>();
/** Submit synchronous state updates; subsequent draws consume the texture on the same GPU queue. */
export function calcMeshColorSmoothingTextureWebGPU(context: WebGPUContext, input: ColorSmoothingInput, options: ColorSmoothingOptions, texture = new WebGPUTextureData()) {
    const { device, width, height, gridDim, samples, bins, gridTransform, type } = prepareSmoothing(context, input, options);
    let pipeline = texturePipelines.get(device);
    if (!pipeline) {
        const module = device.createShaderModule({ label: 'molstar-smoothing-texture', code: textureShader });
        pipeline = device.createComputePipeline({ label: 'molstar-smoothing-texture', layout: 'auto', compute: { module, entryPoint: 'compute' } });
        texturePipelines.set(device, pipeline);
    }
    const output = device.createTexture({ label: 'molstar-smoothed-grid', size: [width, height], format: 'rgba8unorm', usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC });
    const buffers: GPUBuffer[] = [];
    try {
        const params = new Uint32Array(16), x = Math.min(Math.ceil(width * height / 64), device.limits.maxComputeWorkgroupsPerDimension);
        params.set([...gridDim, input.itemSize]); params.set([0, width * height, x, bins.span], 4); params.set(bins.dimensions, 8); params.set([width, height], 12);
        const inputs = [params, samples, bins.ranges, bins.indices].map((a, i) => { const b = context.createBuffer(a, i === 0 ? GPUBufferUsage.UNIFORM : GPUBufferUsage.STORAGE); buffers.push(b); return b; });
        const bindGroup = device.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [...inputs.map((buffer, binding) => ({ binding, resource: { buffer } })), { binding: 4, resource: output.createView() }] });
        const encoder = device.createCommandEncoder({ label: 'molstar-smoothing-texture' });
        const pass = encoder.beginComputePass(); pass.setPipeline(pipeline); pass.setBindGroup(0, bindGroup); pass.dispatchWorkgroups(x, Math.ceil(width * height / (x * 64))); pass.end();
        device.queue.submit([encoder.finish()]); context.stats.computeDispatches++;
        texture.loadGPU(device, output);
    } catch (error) { output.destroy(); throw error; } finally { buffers.forEach(b => b.destroy()); }
    return { kind: 'volume' as const, texture, gridTexDim: Vec2.create(width, height), gridDim, gridTransform, type };
}

export function applyMeshOverlaySmoothingWebGPU(context: WebGPUContext, values: MeshValues | TextureMeshValues, name: 'Overpaint' | 'Transparency' | 'Emissive' | 'Substance', options: ColorSmoothingOptions, texture?: WebGPUTextureData) {
    const props = {
        Overpaint: { type: values.dOverpaintType, colors: values.tOverpaint, grid: values.tOverpaintGrid, texDim: values.uOverpaintTexDim, dim: values.uOverpaintGridDim, transform: values.uOverpaintGridTransform, itemSize: 4 as const },
        Transparency: { type: values.dTransparencyType, colors: values.tTransparency, grid: values.tTransparencyGrid, texDim: values.uTransparencyTexDim, dim: values.uTransparencyGridDim, transform: values.uTransparencyGridTransform, itemSize: 1 as const },
        Emissive: { type: values.dEmissiveType, colors: values.tEmissive, grid: values.tEmissiveGrid, texDim: values.uEmissiveTexDim, dim: values.uEmissiveGridDim, transform: values.uEmissiveGridTransform, itemSize: 1 as const },
        Substance: { type: values.dSubstanceType, colors: values.tSubstance, grid: values.tSubstanceGrid, texDim: values.uSubstanceTexDim, dim: values.uSubstanceGridDim, transform: values.uSubstanceGridTransform, itemSize: 4 as const },
    }[name];
    if (props.type.ref.value !== 'groupInstance') return;
    const data = calcMeshColorSmoothingTextureWebGPU(context, {
        vertexCount: values.uVertexCount.ref.value, instanceCount: values.uInstanceCount.ref.value, groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value, instanceBuffer: values.aInstance.ref.value, ...smoothingGeometry(values),
        colorData: props.colors.ref.value, colorType: 'groupInstance', boundingSphere: values.boundingSphere.ref.value, invariantBoundingSphere: values.invariantBoundingSphere.ref.value, itemSize: props.itemSize,
    }, options, texture);
    ValueCell.updateIfChanged(props.type, 'volumeInstance'); ValueCell.update(props.grid, data.texture);
    ValueCell.update(props.texDim, data.gridTexDim); ValueCell.update(props.dim, data.gridDim); ValueCell.update(props.transform, data.gridTransform);
}
