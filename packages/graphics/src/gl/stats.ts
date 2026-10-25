/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export function createGraphicsStats() {
  return {
    /**
     * Live resources by kind; on WebGPU `attribute` and `elements` are vertex and index buffers, `pixelPack` the staging
     * buffers of reads in flight, `uniform` the uniform buffers (ring chunks, one per arena), `shader` the shader
     * modules (one per program and per kernel pipeline)
     */
    resourceCounts: {
      attribute: 0,
      elements: 0,
      pixelPack: 0,
      uniform: 0,
      framebuffer: 0,
      program: 0,
      renderbuffer: 0,
      shader: 0,
      texture: 0,
      cubeTexture: 0,
      /** WebGL only */
      vertexArray: 0,
      /** WebGPU only, as `kernel` (compute pipelines) */
      storageBuffer: 0,
      kernel: 0,
    },

    /**
     * Draws skipped while their pipeline compiles (WebGPU, in the background), counted since the context was made: a
     * pass that skipped draws compares it before and after
     */
    skippedDraws: 0,

    drawCount: 0,
    instanceCount: 0,
    instancedDrawCount: 0,

    calls: {
      drawInstanced: 0,
      drawInstancedBase: 0,
      multiDrawInstancedBase: 0,
      counts: 0,
    },

    culled: {
      lod: 0,
      frustum: 0,
      occlusion: 0,
    },

    cacheHits: {
      cull: 0,
    },

    /**
     * passes drawn from draw lists (`RendererProps.drawLists`): replayed (`unplanned`: of them, without planning, nothing
     * the plan read changed; `flushed`: with changed data written first), and recorded in part or whole (or drawn,
     * where a recording failed; `chunksRecorded`: the chunks of their objects recorded; `chunksDirect`: drawn without
     * their recordings while changing every frame)
     */
    drawLists: {
      replayed: 0,
      unplanned: 0,
      flushed: 0,
      recorded: 0,
      chunksRecorded: 0,
      chunksDirect: 0,
    },
  };
}

export type GraphicsStats = ReturnType<typeof createGraphicsStats>;

export type ResourceByteCounts = {
  texture: number;
  cubeTexture: number;
  attribute: number;
  elements: number;
  pixelPack: number;
  uniform: number;
  renderbuffer: number;
  /** storage buffers of compute kernels and render items (WebGPU) */
  storageBuffer: number;
};
