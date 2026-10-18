/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * The state set on a render pass encoder, so that draws skip calls that repeat it (consecutive draws of one pipeline,
 * the global bind group of a pass, the viewport). Commands recorded into the encoder directly must `invalidate` it.
 */
export interface PassState {
  readonly encoder: GPURenderPassEncoder;
  setPipeline: (pipeline: GPURenderPipeline) => void;
  /** `offsets` are compared by value (callers may reuse the array) */
  setBindGroup: (index: number, group: GPUBindGroup, offsets?: readonly number[]) => void;
  setVertexBuffer: (slot: number, buffer: GPUBuffer) => void;
  setIndexBuffer: (buffer: GPUBuffer, format: GPUIndexFormat) => void;
  setViewport: (x: number, y: number, width: number, height: number) => void;
  setScissorRect: (x: number, y: number, width: number, height: number) => void;
  setStencilReference: (reference: number) => void;
  /** forgets the state, after commands recorded into `encoder` directly */
  invalidate: () => void;
}

const NoOffsets: readonly number[] = [];

export function createPassState(encoder: GPURenderPassEncoder): PassState {
  let pipeline: GPURenderPipeline | undefined;
  const groups: (GPUBindGroup | undefined)[] = [];
  const groupOffsets: number[][] = [];
  const vertexBuffers: (GPUBuffer | undefined)[] = [];
  let indexBuffer: GPUBuffer | undefined;
  const viewport = [NaN, NaN, NaN, NaN];
  const scissor = [NaN, NaN, NaN, NaN];
  let stencilReference = NaN;

  const sameRect = (rect: number[], x: number, y: number, width: number, height: number) =>
    rect[0] === x && rect[1] === y && rect[2] === width && rect[3] === height;
  const setRect = (rect: number[], x: number, y: number, width: number, height: number) => {
    rect[0] = x;
    rect[1] = y;
    rect[2] = width;
    rect[3] = height;
  };

  return {
    encoder,
    setPipeline: (p) => {
      if (p === pipeline) return;
      pipeline = p;
      encoder.setPipeline(p);
    },
    setBindGroup: (index, group, offsets = NoOffsets) => {
      let last = groupOffsets[index];
      if (groups[index] === group && last && last.length === offsets.length) {
        let same = true;
        for (let i = 0, il = offsets.length; i < il; ++i) {
          if (last[i] !== offsets[i]) {
            same = false;
            break;
          }
        }
        if (same) return;
      }
      groups[index] = group;
      if (!last) groupOffsets[index] = last = [];
      last.length = offsets.length;
      for (let i = 0, il = offsets.length; i < il; ++i) last[i] = offsets[i];
      encoder.setBindGroup(index, group, offsets);
    },
    setVertexBuffer: (slot, buffer) => {
      if (vertexBuffers[slot] === buffer) return;
      vertexBuffers[slot] = buffer;
      encoder.setVertexBuffer(slot, buffer);
    },
    setIndexBuffer: (buffer, format) => {
      if (indexBuffer === buffer) return;
      indexBuffer = buffer;
      encoder.setIndexBuffer(buffer, format);
    },
    setViewport: (x, y, width, height) => {
      if (sameRect(viewport, x, y, width, height)) return;
      setRect(viewport, x, y, width, height);
      encoder.setViewport(x, y, width, height, 0, 1);
    },
    setScissorRect: (x, y, width, height) => {
      if (sameRect(scissor, x, y, width, height)) return;
      setRect(scissor, x, y, width, height);
      encoder.setScissorRect(x, y, width, height);
    },
    setStencilReference: (reference) => {
      if (reference === stencilReference) return;
      stencilReference = reference;
      encoder.setStencilReference(reference);
    },
    invalidate: () => {
      pipeline = undefined;
      groups.length = 0;
      groupOffsets.length = 0;
      vertexBuffers.length = 0;
      indexBuffer = undefined;
      setRect(viewport, NaN, NaN, NaN, NaN);
      setRect(scissor, NaN, NaN, NaN, NaN);
      stencilReference = NaN;
    },
  };
}
