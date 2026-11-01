/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BufferUsage, MapMode, ShaderStage, TextureUsage } from './constants.js';

/** the depth of the probe's triangles, inside the viewport's depth range */
const ProbeDepth = 0.75;

const ProbeShader = /* wgsl */ `
@group(0) @binding(0) var<uniform> uA: vec4<f32>;
@group(1) @binding(0) var<uniform> uB: vec4<f32>;

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
};

@vertex
fn vertexMain(@builtin(vertex_index) i: u32) -> VertexOutput {
    let p = array(vec2<f32>(-1.0, -1.0), vec2<f32>(3.0, -1.0), vec2<f32>(-1.0, 3.0));
    return VertexOutput(vec4<f32>(p[i], ${ProbeDepth}, 1.0));
}

@fragment
fn fragmentPlain() -> @location(0) vec4<f32> {
    return vec4<f32>(0.5) * uA * uB;
}

struct DepthOutput {
    @location(0) color: vec4<f32>,
    @builtin(frag_depth) fragDepth: f32,
};

@fragment
fn fragmentDepth(in: VertexOutput) -> DepthOutput {
    return DepthOutput(uA * uB, in.position.z);
}
`;

/**
 * Whether a draw writing `frag_depth` keeps its depth when replayed from a render bundle. Safari 26 (iOS) loses it
 * (the draw fails its depth test, invisible) where the draw follows another pipeline in the same bundle and the
 * pipeline layout has two or more uniform buffers at dynamic offsets, as our programs (globals, material, object) and
 * draw lists (several programs per bundle) do; WebKit clamps the written depth to the viewport's depth range, passed to
 * the shader (WebKit bug 264449). One draw of each pipeline into a 1x1 target, the color read back: the plain one's
 * half, replaced at the same depth (`greater-equal`) where the depth is kept; false where the probe fails.
 */
export async function testBundleFragDepth(device: GPUDevice): Promise<boolean> {
  const destroyable: { destroy(): void }[] = [];
  try {
    // the uniforms (1) at offset 256, bound there, zeros below
    const uniforms = device.createBuffer({ size: 512, usage: BufferUsage.UNIFORM | BufferUsage.COPY_DST });
    destroyable.push(uniforms);
    device.queue.writeBuffer(uniforms, 256, new Float32Array([1, 1, 1, 1]));
    // as reproduced on the device: one uniform buffer at a dynamic offset in each of two groups
    const groupLayout = device.createBindGroupLayout({
      entries: [
        {
          binding: 0,
          visibility: ShaderStage.FRAGMENT,
          buffer: { type: 'uniform', hasDynamicOffset: true, minBindingSize: 16 },
        },
      ],
    });
    const group = device.createBindGroup({
      layout: groupLayout,
      entries: [{ binding: 0, resource: { buffer: uniforms, size: 16 } }],
    });
    const layout = device.createPipelineLayout({ bindGroupLayouts: [groupLayout, groupLayout] });
    const module = device.createShaderModule({ code: ProbeShader, label: 'bundle frag depth probe' });
    const colorFormat: GPUTextureFormat = 'rgba8unorm';
    const depthFormat: GPUTextureFormat = 'depth24plus-stencil8';
    const pipeline = (entryPoint: string, depthCompare: GPUCompareFunction) =>
      device.createRenderPipelineAsync({
        layout,
        vertex: { module, entryPoint: 'vertexMain' },
        fragment: { module, entryPoint, targets: [{ format: colorFormat }] },
        primitive: { topology: 'triangle-list' },
        depthStencil: { format: depthFormat, depthWriteEnabled: true, depthCompare },
        label: 'bundle frag depth probe',
      });
    const [plain, depth] = await Promise.all([
      pipeline('fragmentPlain', 'greater'),
      pipeline('fragmentDepth', 'greater-equal'),
    ]);

    const bundleEncoder = device.createRenderBundleEncoder({
      colorFormats: [colorFormat],
      depthStencilFormat: depthFormat,
    });
    for (const p of [plain, depth]) {
      bundleEncoder.setPipeline(p);
      bundleEncoder.setBindGroup(0, group, [256]);
      bundleEncoder.setBindGroup(1, group, [256]);
      bundleEncoder.draw(3);
    }
    const bundle = bundleEncoder.finish();

    const usage = TextureUsage.RENDER_ATTACHMENT | TextureUsage.COPY_SRC;
    const color = device.createTexture({ size: [1, 1], format: colorFormat, usage });
    const depthTexture = device.createTexture({ size: [1, 1], format: depthFormat, usage });
    // a row of a copy is aligned to 256 bytes
    const readback = device.createBuffer({ size: 256, usage: BufferUsage.COPY_DST | BufferUsage.MAP_READ });
    destroyable.push(color, depthTexture, readback);
    const encoder = device.createCommandEncoder();
    const pass = encoder.beginRenderPass({
      colorAttachments: [{ view: color.createView(), clearValue: [0, 0, 0, 0], loadOp: 'clear', storeOp: 'store' }],
      depthStencilAttachment: {
        view: depthTexture.createView(),
        depthClearValue: 0,
        depthLoadOp: 'clear',
        depthStoreOp: 'store',
        stencilClearValue: 0,
        stencilLoadOp: 'clear',
        stencilStoreOp: 'store',
      },
    });
    pass.setViewport(0, 0, 1, 1, 0, 1);
    pass.executeBundles([bundle]);
    pass.end();
    encoder.copyTextureToBuffer({ texture: color }, { buffer: readback, bytesPerRow: 256 }, [1, 1]);
    device.queue.submit([encoder.finish()]);

    await readback.mapAsync(MapMode.READ);
    const red = new Uint8Array(readback.getMappedRange(), 0, 1)[0];
    readback.unmap();
    return red === 255;
  } catch {
    return false;
  } finally {
    for (const d of destroyable) d.destroy();
  }
}
