/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { WebGPUContext } from '../../../mol-gl/webgpu/context';
import { GPUBufferUsage, GPUMapMode } from '../../../mol-gl/webgpu/compat';
import { RuntimeContext } from '../../../mol-task';
import { AlphaOrbital, CubeGridInfo } from '../data-model';
import { createTextureData, getNormalizedAlpha } from './data';
import { orbitalShader } from './shader.wgsl';

const pipelines = new WeakMap<GPUDevice, Promise<GPUComputePipeline>>();
function getPipeline(device: GPUDevice) {
    let pipeline = pipelines.get(device);
    if (!pipeline) {
        pipeline = (async () => {
            const module = device.createShaderModule({ label: 'molstar-orbitals', code: orbitalShader });
            const errors = (await module.getCompilationInfo()).messages.filter(m => m.type === 'error');
            if (errors.length) throw new Error(`WebGPU orbital shader: ${errors.map(e => e.message).join('; ')}`);
            return device.createComputePipelineAsync({ layout: 'auto', compute: { module, entryPoint: 'compute' } });
        })();
        pipelines.set(device, pipeline);
    }
    return pipeline;
}

/** Spherical orbital/density values in the CPU grid's z-fastest storage order. */
export async function computeOrbitalGridWebGPU(ctx: RuntimeContext, context: WebGPUContext, grid: CubeGridInfo, orbitals: AlphaOrbital[], density = false, options: { maxVoxelsPerBatch?: number } = {}): Promise<Float32Array> {
    if (!density && orbitals.length !== 1) throw new Error('Orbital computation requires exactly one orbital.');
    if (!grid.dimensions.every(d => Number.isSafeInteger(d) && d > 0) || !Number.isSafeInteger(grid.npoints) || grid.npoints !== grid.dimensions[0] * grid.dimensions[1] * grid.dimensions[2] || grid.npoints > 0xffffffff) throw new Error('WebGPU orbital grid exceeds the addressable voxel range.');
    if (![...grid.box.min, ...grid.delta].every(v => Number.isFinite(Math.fround(v)))) throw new Error('WebGPU orbital grid coordinates must be finite f32 values.');
    for (const atom of grid.params.basis.atoms) {
        if (!atom.center.every(v => Number.isFinite(Math.fround(v)))) throw new Error('WebGPU orbital centers must be finite f32 values.');
        for (const shell of atom.shells) {
            if (!shell.exponents.length || !shell.exponents.every(e => Number.isFinite(Math.fround(e)) && Math.fround(e) > 0)) throw new Error('WebGPU orbital exponents must be positive finite f32 values.');
            if (!shell.angularMomentum.every(l => Number.isInteger(l) && l >= 0 && l <= 4) || shell.coefficients.length !== shell.angularMomentum.length || !shell.coefficients.every(c => c.length === shell.exponents.length && c.every(v => Number.isFinite(Math.fround(v))))) throw new Error('Invalid WebGPU orbital shell coefficients or angular momentum.');
        }
    }
    if (!Number.isFinite(grid.params.cutoffThreshold) || grid.params.cutoffThreshold > 1) throw new Error('Invalid WebGPU orbital cutoff threshold.');
    const alphaCount = grid.params.basis.atoms.reduce((n, a) => n + a.shells.reduce((s, shell) => s + shell.angularMomentum.reduce((t, l) => t + 2 * l + 1, 0), 0), 0);
    for (const orbital of orbitals) {
        if (orbital.alpha.length !== alphaCount || !orbital.alpha.every(v => Number.isFinite(Math.fround(v))) || !Number.isFinite(Math.fround(orbital.occupancy))) throw new Error('Invalid WebGPU orbital alpha coefficients or occupancy.');
    }
    const result = new Float32Array(grid.npoints);
    const active = density ? orbitals.filter(o => o.occupancy !== 0) : orbitals;
    if (!active.length || !alphaCount) return result;
    const data = createTextureData(grid, active[0]);
    const alpha = new Float32Array(alphaCount * active.length);
    active.forEach((o, i) => alpha.set(getNormalizedAlpha(grid.params.basis, o.alpha, grid.params.sphericalOrder), i * alphaCount));
    const inputs = [data.tCenters.array, data.tInfo.array, data.tCoeff.array, alpha, Float32Array.from(active.map(o => o.occupancy))];
    const { device } = context;
    const limit = Math.min(device.limits.maxStorageBufferBindingSize, device.limits.maxBufferSize);
    if (inputs.some(a => a.byteLength > limit || !a.every(Number.isFinite))) throw new Error('WebGPU orbital basis exceeds device limits or finite f32 range.');
    const requested = options.maxVoxelsPerBatch ?? 262144;
    if (!Number.isSafeInteger(requested) || requested < 1) throw new Error('WebGPU orbital batch size must be a positive integer.');
    const batchSize = Math.min(grid.npoints, requested, Math.floor(limit / 4), device.limits.maxComputeWorkgroupsPerDimension ** 2 * 64);
    const pipeline = await getPipeline(device);
    const buffers: GPUBuffer[] = [];
    try {
        const bindings = inputs.map(a => { const b = context.createBuffer(a, GPUBufferUsage.STORAGE); buffers.push(b); return b; });
        const uniform = device.createBuffer({ size: 80, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST }); buffers.push(uniform);
        const output = device.createBuffer({ size: batchSize * 4, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC }); buffers.push(output);
        const readback = device.createBuffer({ size: batchSize * 4, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST }); buffers.push(readback);
        const bindGroup = device.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [uniform, ...bindings, output].map((buffer, binding) => ({ binding, resource: { buffer } })) });
        const parameters = new Float32Array(20), uints = new Uint32Array(parameters.buffer);
        uints.set(grid.dimensions); parameters.set(grid.box.min, 4); parameters.set(grid.delta, 8);
        uints.set([data.uNCenters, active.length, density ? 1 : 0, alphaCount], 16);
        for (let base = 0; base < grid.npoints; base += batchSize) {
            const count = Math.min(batchSize, grid.npoints - base);
            const x = Math.min(Math.ceil(count / 64), device.limits.maxComputeWorkgroupsPerDimension);
            uints.set([base, count, x, 0], 12); device.queue.writeBuffer(uniform, 0, parameters);
            const encoder = device.createCommandEncoder({ label: 'molstar-orbital-grid' });
            const pass = encoder.beginComputePass(); pass.setPipeline(pipeline); pass.setBindGroup(0, bindGroup);
            pass.dispatchWorkgroups(x, Math.ceil(count / (x * 64))); pass.end();
            encoder.copyBufferToBuffer(output, 0, readback, 0, count * 4); device.queue.submit([encoder.finish()]); context.stats.computeDispatches++;
            await readback.mapAsync(GPUMapMode.READ, 0, count * 4);
            result.set(new Float32Array(readback.getMappedRange(0, count * 4)), base); readback.unmap();
            if (ctx.shouldUpdate) await ctx.update({ message: 'Computing orbital grid on WebGPU', current: base + count, max: grid.npoints });
        }
    } finally { buffers.forEach(b => b.destroy()); }
    return result;
}
