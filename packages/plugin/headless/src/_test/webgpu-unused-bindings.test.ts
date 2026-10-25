/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Canvas3DParams, type Canvas3DPrecompileProps } from '@molstar/graphics/canvas3d/canvas3d';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `${Array.from({ length: 40 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const residue = String(i + 1).padStart(4);
  return `ATOM  ${serial}  CA  ALA A${residue}    ${f(Math.cos(i) * 8)}${f(Math.sin(i) * 8)}${f(i * 1.5)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

type Layout = { entries: GPUBindGroupLayoutEntry[] };
type Recorded = { label: string; code: string; descriptor: GPURenderPipelineDescriptor; groups: Layout[] };

/** The props of the frames whose pipelines are checked: the defaults and each pass turned on */
function getPropsVariants(): Canvas3DPrecompileProps[] {
  const defaults = PD.getDefaultValues(Canvas3DParams);
  const pp = Canvas3DParams.postprocessing.params as Record<string, any>;
  const on = (name: string, variant = 'on') => ({
    name: variant,
    params: pp[name].map(variant).defaultValue,
  });
  const postprocessing = (changes: Record<string, unknown>) => ({
    ...defaults,
    postprocessing: { ...defaults.postprocessing, ...changes },
  });
  return [
    { ...defaults, transparency: 'wboit' },
    { ...defaults, transparency: 'dpoit' },
    { ...defaults, transparency: 'blended' },
    postprocessing({ occlusion: on('occlusion') }),
    postprocessing({ shadow: on('shadow') }),
    postprocessing({ outline: on('outline') }),
    postprocessing({ dof: on('dof') }),
    postprocessing({ bloom: { name: 'off', params: {} } }),
    postprocessing({ antialiasing: on('antialiasing', 'smaa') }),
    postprocessing({ sharpening: on('sharpening') }),
    postprocessing({ enabled: false }),
    { ...defaults, multiSample: { ...defaults.multiSample, mode: 'on' } },
    { ...defaults, illumination: { ...defaults.illumination, enabled: true } },
  ] as Canvas3DPrecompileProps[];
}

// a program's bind group layouts bind only what its shader statically uses: what a variant turns off is left out of
// its source (`WgslProgram.source`), so that unused targets are neither bound nor allocated; WebGPU is the judge:
// without a binding the shader doesn't use, the pipeline still validates
describe.skipIf(!WebGPU.module)('WebGPU bindings', () => {
  it('bind only what the shaders use', async () => {
    const g = WebGPU.module!.globals as any;
    const codes = new WeakMap<GPUShaderModule, string>();
    const groupLayouts = new WeakMap<GPUBindGroupLayout, Layout>();
    const pipelineLayouts = new WeakMap<GPUPipelineLayout, Layout[]>();
    const recorded = new Map<GPUShaderModule, Recorded>();
    const proto = g.GPUDevice.prototype;
    const original = {
      createShaderModule: proto.createShaderModule,
      createBindGroupLayout: proto.createBindGroupLayout,
      createPipelineLayout: proto.createPipelineLayout,
      createRenderPipeline: proto.createRenderPipeline,
      createRenderPipelineAsync: proto.createRenderPipelineAsync,
    };
    const record = (d: GPURenderPipelineDescriptor) => {
      const groups = d.layout !== 'auto' ? pipelineLayouts.get(d.layout) : undefined;
      const code = codes.get(d.vertex.module);
      if (!groups || !code || !d.label) return;
      // one per program variant, a shader module each (pipelines of other states and targets bind the same)
      const label = d.label.split(' ').slice(0, 2).join(' ');
      if (!recorded.has(d.vertex.module)) recorded.set(d.vertex.module, { label, code, descriptor: d, groups });
    };
    proto.createShaderModule = function (this: GPUDevice, d: GPUShaderModuleDescriptor) {
      const m = original.createShaderModule.call(this, d);
      codes.set(m, d.code);
      return m;
    };
    proto.createBindGroupLayout = function (this: GPUDevice, d: GPUBindGroupLayoutDescriptor) {
      const l = original.createBindGroupLayout.call(this, d);
      groupLayouts.set(l, { entries: [...d.entries] });
      return l;
    };
    proto.createPipelineLayout = function (this: GPUDevice, d: GPUPipelineLayoutDescriptor) {
      const l = original.createPipelineLayout.call(this, d);
      pipelineLayouts.set(
        l,
        [...d.bindGroupLayouts].map((b) => groupLayouts.get(b as GPUBindGroupLayout)!),
      );
      return l;
    };
    proto.createRenderPipeline = function (this: GPUDevice, d: GPURenderPipelineDescriptor) {
      record(d);
      return original.createRenderPipeline.call(this, d);
    };
    proto.createRenderPipelineAsync = function (this: GPUDevice, d: GPURenderPipelineDescriptor) {
      record(d);
      return original.createRenderPipelineAsync.call(this, d);
    };

    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 64, height: 64 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const { representation } = plugin.builders.structure;
      await representation.addRepresentation(structure, { type: 'spacefill', typeParams: { alpha: 0.5 } });
      await representation.addRepresentation(structure, { type: 'cartoon' });
      const canvas = plugin.canvas3d!;
      canvas.commit(true);
      await plugin.getImageRaw({ width: 64, height: 64 });
      await canvas.precompile(canvas.getPrecompileList(), getPropsVariants());
      await canvas.prepare({ interaction: true });
      await canvas.graphics.resources.pending.whenReady();

      const device = (canvas.graphics as unknown as { device: GPUDevice }).device;
      const unused: string[] = [];
      for (const { label, code, descriptor, groups } of recorded.values()) {
        for (let i = 0; i < groups.length; ++i) {
          for (const entry of groups[i].entries) {
            const without = groups.map((group, j) =>
              original.createBindGroupLayout.call(device, {
                entries: j === i ? group.entries.filter((e) => e !== entry) : group.entries,
              }),
            );
            const layout = original.createPipelineLayout.call(device, { bindGroupLayouts: without });
            device.pushErrorScope('validation');
            original.createRenderPipeline.call(device, { ...descriptor, layout });
            const error = await device.popErrorScope();
            if (error) continue;
            const declaration = new RegExp(
              `@group\\(${i}\\)\\s*@binding\\(${entry.binding}\\)\\s*var(?:<[^>]*>)?\\s+(\\w+)`,
            );
            const defines = [...code.matchAll(/const (d\w+)\s*=\s*(\w+)/g)].map((m) => `${m[1]}=${m[2]}`).join(' ');
            unused.push(
              `${label}: ${declaration.exec(code)?.[1] ?? `group ${i} binding ${entry.binding}`} (${defines})`,
            );
          }
        }
      }
      expect(recorded.size).toBeGreaterThan(20);
      expect(unused).toEqual([]);
    } finally {
      Object.assign(proto, original);
      plugin.dispose();
    }
  });
});
