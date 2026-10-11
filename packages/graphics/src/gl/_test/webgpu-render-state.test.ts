/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BlendState, type RenderState } from '../render-state.js';
import { getRenderStateId } from '../webgpu/render-state.js';

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
        expect(getRenderStateId({ stencil })).not.toBe(
            getRenderStateId({ stencil: { ...stencil, writeMask: 0 } }),
        );
    });
});
