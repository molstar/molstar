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
  /** the render pass, or a render bundle being recorded (`createBundlePassState`) */
  readonly encoder: GPURenderPassEncoder | GPURenderBundleEncoder;
  setPipeline: (pipeline: GPURenderPipeline) => void;
  /** `offsets` are compared by value (callers may reuse the array) */
  setBindGroup: (index: number, group: GPUBindGroup, offsets?: readonly number[]) => void;
  /** binds `size` bytes of `buffer` from `offset` (ranges of shared buffers, `BufferPages`) */
  setVertexBuffer: (slot: number, buffer: GPUBuffer, offset: number, size: number) => void;
  setIndexBuffer: (buffer: GPUBuffer, format: GPUIndexFormat, offset: number, size: number) => void;
  setViewport: (x: number, y: number, width: number, height: number) => void;
  setScissorRect: (x: number, y: number, width: number, height: number) => void;
  setStencilReference: (reference: number) => void;
  /** forgets the state, after commands recorded into `encoder` directly */
  invalidate: () => void;
}

const NoOffsets: readonly number[] = [];

/**
 * An encoder that only lists what the commands bind, as text (`GraphicsContext.captureDraws`): pipelines, bind groups
 * with offsets, buffers and draw arguments, objects by an id per object; viewport, scissor and stencil reference are the
 * pass's (not of a recording) and left out
 */
export function createCaptureEncoder(tokens: string[], getId: (o: object) => number): GPURenderPassEncoder {
  const encoder = {
    setPipeline: (p: GPURenderPipeline) => tokens.push(`p${getId(p)}`),
    setBindGroup: (index: number, group: GPUBindGroup | null, offsets?: Iterable<number>) =>
      tokens.push(`g${index}:${group ? getId(group) : -1}:${offsets ? Array.from(offsets).join(',') : ''}`),
    setVertexBuffer: (slot: number, buffer: GPUBuffer | null, offset = 0, size?: number) =>
      tokens.push(`v${slot}:${buffer ? getId(buffer) : -1}:${offset}:${size ?? ''}`),
    setIndexBuffer: (buffer: GPUBuffer, format: GPUIndexFormat, offset = 0, size?: number) =>
      tokens.push(`i${getId(buffer)}:${format}:${offset}:${size ?? ''}`),
    draw: (...args: number[]) => tokens.push(`d${args.join(',')}`),
    drawIndexed: (...args: number[]) => tokens.push(`x${args.join(',')}`),
    drawIndirect: (buffer: GPUBuffer, offset: number) => tokens.push(`r${getId(buffer)}:${offset}`),
    drawIndexedIndirect: (buffer: GPUBuffer, offset: number) => tokens.push(`s${getId(buffer)}:${offset}`),
    setViewport: () => {},
    setScissorRect: () => {},
    setStencilReference: () => {},
    setBlendConstant: () => {},
  };
  return encoder as unknown as GPURenderPassEncoder;
}

/**
 * The state of a render bundle being recorded for the pass of `pass`: pipelines, bind groups and buffers go into the
 * bundle; viewport, scissor and stencil reference are the pass's (a bundle has none, `executeBundles` keeps the pass's),
 * set on `pass` once; a different value later in the recording calls `onConflict` (the bundle can't replay it)
 */
export function createBundlePassState(
  encoder: GPURenderBundleEncoder,
  pass: PassState,
  onConflict: (what: string) => void,
): PassState {
  const state = createPassState(encoder);
  const viewport = [NaN, NaN, NaN, NaN];
  const scissor = [NaN, NaN, NaN, NaN];
  let stencilReference = NaN;
  const fixed = (rect: number[], x: number, y: number, width: number, height: number, what: string) => {
    if (Number.isNaN(rect[0])) {
      rect[0] = x;
      rect[1] = y;
      rect[2] = width;
      rect[3] = height;
      return true;
    }
    if (rect[0] !== x || rect[1] !== y || rect[2] !== width || rect[3] !== height) onConflict(what);
    return false;
  };
  return {
    ...state,
    encoder,
    setViewport: (x, y, width, height) => {
      if (fixed(viewport, x, y, width, height, 'viewport')) pass.setViewport(x, y, width, height);
    },
    setScissorRect: (x, y, width, height) => {
      if (fixed(scissor, x, y, width, height, 'scissor')) pass.setScissorRect(x, y, width, height);
    },
    setStencilReference: (reference) => {
      if (Number.isNaN(stencilReference)) {
        stencilReference = reference;
        pass.setStencilReference(reference);
      } else if (reference !== stencilReference) onConflict('stencil reference');
    },
  };
}

/** The state of a pass, or of a bundle (whose viewport, scissor and stencil `createBundlePassState` replaces) */
export function createPassState(encoder: GPURenderPassEncoder | GPURenderBundleEncoder): PassState {
  let pipeline: GPURenderPipeline | undefined;
  const groups: (GPUBindGroup | undefined)[] = [];
  const groupOffsets: number[][] = [];
  const vertexBuffers: (GPUBuffer | undefined)[] = [];
  const vertexOffsets: number[] = [];
  let indexBuffer: GPUBuffer | undefined;
  let indexOffset = NaN;
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
    setVertexBuffer: (slot, buffer, offset, size) => {
      // a range is bound with its size: the same buffer and offset bind the same range
      if (vertexBuffers[slot] === buffer && vertexOffsets[slot] === offset) return;
      vertexBuffers[slot] = buffer;
      vertexOffsets[slot] = offset;
      encoder.setVertexBuffer(slot, buffer, offset, size);
    },
    setIndexBuffer: (buffer, format, offset, size) => {
      if (indexBuffer === buffer && indexOffset === offset) return;
      indexBuffer = buffer;
      indexOffset = offset;
      encoder.setIndexBuffer(buffer, format, offset, size);
    },
    setViewport: (x, y, width, height) => {
      if (sameRect(viewport, x, y, width, height)) return;
      setRect(viewport, x, y, width, height);
      (encoder as GPURenderPassEncoder).setViewport(x, y, width, height, 0, 1);
    },
    setScissorRect: (x, y, width, height) => {
      if (sameRect(scissor, x, y, width, height)) return;
      setRect(scissor, x, y, width, height);
      (encoder as GPURenderPassEncoder).setScissorRect(x, y, width, height);
    },
    setStencilReference: (reference) => {
      if (reference === stencilReference) return;
      stencilReference = reference;
      (encoder as GPURenderPassEncoder).setStencilReference(reference);
    },
    invalidate: () => {
      pipeline = undefined;
      groups.length = 0;
      groupOffsets.length = 0;
      vertexBuffers.length = 0;
      vertexOffsets.length = 0;
      indexBuffer = undefined;
      indexOffset = NaN;
      setRect(viewport, NaN, NaN, NaN, NaN);
      setRect(scissor, NaN, NaN, NaN, NaN);
      stencilReference = NaN;
    },
  };
}
