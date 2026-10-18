/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { BlendState, intersectTargetRect, type RenderState } from '../render-state.js';
import { getDepthStencilState, getRenderStateId } from '../webgpu/render-state.js';

describe('webgpu render state ids', () => {
  const stencil = {
    front: { compare: 'always', failOp: 'keep', depthFailOp: 'keep', passOp: 'increment-wrap' },
    back: { compare: 'always', failOp: 'keep', depthFailOp: 'keep', passOp: 'decrement-wrap' },
    readMask: 0xff,
    writeMask: 0xff,
    reference: 0,
  } as const;

  it('are equal for states with equal content, also built per draw', () => {
    const base: RenderState = { depth: { compare: 'less', write: true }, cull: 'back' };
    const a: RenderState = { ...base, colorWrite: [false, false, false, false], stencil };
    const b: RenderState = { stencil: { ...stencil }, colorWrite: [false, false, false, false], ...base };
    expect(getRenderStateId(a)).toBe(getRenderStateId(b));
  });

  it('ignore the stencil reference (dynamic state)', () => {
    expect(getRenderStateId({ stencil })).toBe(getRenderStateId({ stencil: { ...stencil, reference: 1 } }));
  });

  it('differ for different content', () => {
    const blend = BlendState.Max;
    expect(getRenderStateId({ blend })).not.toBe(getRenderStateId({}));
    expect(getRenderStateId({ cull: 'back' })).not.toBe(getRenderStateId({ cull: 'front' }));
    expect(getRenderStateId({ stencil })).not.toBe(getRenderStateId({ stencil: { ...stencil, writeMask: 0 } }));
  });
});

describe('target rects', () => {
  const rect = (x: number, y: number, width: number, height: number) =>
    intersectTargetRect({ x: 0, y: 0, width: 0, height: 0 }, x, y, width, height, 100, 100);

  it('keep rects inside the target', () => {
    expect(rect(10, 20, 30, 40)).toEqual({ x: 10, y: 20, width: 30, height: 40 });
  });

  it('are the part inside the target, on every side', () => {
    expect(rect(-10, -10, 50, 50)).toEqual({ x: 0, y: 0, width: 40, height: 40 });
    expect(rect(80, 80, 50, 50)).toEqual({ x: 80, y: 80, width: 20, height: 20 });
  });

  it('are empty outside the target', () => {
    expect(rect(120, 0, 10, 10)).toEqual({ x: 100, y: 0, width: 0, height: 10 });
    expect(rect(0, -30, 10, 10)).toEqual({ x: 0, y: 0, width: 10, height: 0 });
  });
});

describe('depth stencil states', () => {
  const state: RenderState = {
    depth: { compare: 'greater', write: true },
    stencil: {
      front: { compare: 'always', failOp: 'keep', depthFailOp: 'keep', passOp: 'replace' },
      back: { compare: 'always', failOp: 'keep', depthFailOp: 'keep', passOp: 'replace' },
      readMask: 0xff,
      writeMask: 0xff,
      reference: 1,
    },
  };

  it('write as the render state says', () => {
    const s = getDepthStencilState(state, 'depth24plus-stencil8', true);
    expect(s.depthWriteEnabled).toBe(true);
    expect(s.stencilWriteMask).toBe(0xff);
  });

  it('mask depth and stencil writes with a read-only depth attachment, the tests stay', () => {
    const s = getDepthStencilState(state, 'depth24plus-stencil8', true, true);
    expect(s.depthWriteEnabled).toBe(false);
    expect(s.stencilWriteMask).toBe(0);
    expect(s.depthCompare).toBe('greater');
    expect(s.stencilFront?.compare).toBe('always');
  });
});
