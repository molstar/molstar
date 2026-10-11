/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** Mip level count of a full chain for a `width` by `height` texture */
export function getMipLevelCount(width: number, height: number) {
  return Math.floor(Math.log2(Math.max(1, width, height))) + 1;
}

const MipmapShader = /* wgsl */ `
struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) uv: vec2<f32>,
};

@vertex
fn vertexMain(@builtin(vertex_index) i: u32) -> VertexOutput {
    // one triangle covering the target
    let p = array<vec2<f32>, 3>(vec2<f32>(-1.0, -1.0), vec2<f32>(3.0, -1.0), vec2<f32>(-1.0, 3.0));
    var out: VertexOutput;
    out.position = vec4<f32>(p[i], 0.0, 1.0);
    out.uv = p[i] * vec2<f32>(0.5, -0.5) + 0.5;
    return out;
}

@group(0) @binding(0) var sLinear: sampler;
@group(0) @binding(1) var tSource: texture_2d<f32>;

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    return textureSampleLevel(tSource, sLinear, in.uv, 0.0);
}
`;

/**
 * Fills the mip levels below level 0 of 2d textures (all array layers, e.g. the faces of cube textures) by linearly
 * downsampling the level above, as `generateMipmap` on WebGL; the textures need render attachment usage
 */
export function createMipmapGenerator(device: GPUDevice) {
  let module: GPUShaderModule | undefined;
  let sampler: GPUSampler | undefined;
  const pipelines = new Map<GPUTextureFormat, GPURenderPipeline>();

  function getPipeline(format: GPUTextureFormat) {
    let pipeline = pipelines.get(format);
    if (!pipeline) {
      module ??= device.createShaderModule({ code: MipmapShader, label: 'mipmap' });
      pipeline = device.createRenderPipeline({
        layout: 'auto',
        vertex: { module, entryPoint: 'vertexMain' },
        fragment: { module, entryPoint: 'fragmentMain', targets: [{ format }] },
        primitive: { topology: 'triangle-list' },
        label: `mipmap ${format}`,
      });
      pipelines.set(format, pipeline);
    }
    return pipeline;
  }

  return (encoder: GPUCommandEncoder, texture: GPUTexture) => {
    if (texture.mipLevelCount < 2) return;
    sampler ??= device.createSampler({ magFilter: 'linear', minFilter: 'linear', label: 'mipmap' });
    const pipeline = getPipeline(texture.format);
    for (let layer = 0; layer < texture.depthOrArrayLayers; ++layer) {
      for (let level = 1; level < texture.mipLevelCount; ++level) {
        const view = (baseMipLevel: number) =>
          texture.createView({
            dimension: '2d',
            baseMipLevel,
            mipLevelCount: 1,
            baseArrayLayer: layer,
            arrayLayerCount: 1,
          });
        const bindGroup = device.createBindGroup({
          layout: pipeline.getBindGroupLayout(0),
          entries: [
            { binding: 0, resource: sampler },
            { binding: 1, resource: view(level - 1) },
          ],
        });
        const pass = encoder.beginRenderPass({
          colorAttachments: [{ view: view(level), loadOp: 'clear', storeOp: 'store', clearValue: [0, 0, 0, 0] }],
          label: 'mipmap',
        });
        pass.setPipeline(pipeline);
        pass.setBindGroup(0, bindGroup);
        pass.draw(3);
        pass.end();
      }
    }
  };
}

export type MipmapGenerator = ReturnType<typeof createMipmapGenerator>;
