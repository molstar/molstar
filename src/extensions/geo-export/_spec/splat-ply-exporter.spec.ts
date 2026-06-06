/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 */

import { Color } from '../../../mol-util/color/color';
import { encodeSplat } from '../splat-ply-exporter';

const C0 = 0.28209479177387814;

describe('splat-ply encodeSplat (INRIA 3DGS conventions)', () => {
    it('stores scale as log(radius)', () => {
        expect(encodeSplat(2, 1, Color(0xffffff)).logScale).toBeCloseTo(Math.log(2), 10);
        expect(encodeSplat(1, 1, Color(0xffffff)).logScale).toBeCloseTo(0, 10);
    });

    it('stores opacity as logit(alpha), clamped', () => {
        expect(encodeSplat(1, 0.5, Color(0)).logitOpacity).toBeCloseTo(0, 10); // logit(0.5) = 0
        expect(encodeSplat(1, 1, Color(0)).logitOpacity).toBeGreaterThan(0); // fully opaque -> positive (clamped)
        expect(encodeSplat(1, 0, Color(0)).logitOpacity).toBeLessThan(0); // fully transparent -> negative (clamped)
    });

    it('stores color as SH DC term f_dc = (c - 0.5) / C0', () => {
        for (const v of encodeSplat(1, 1, Color(0xffffff)).fDc) expect(v).toBeCloseTo(0.5 / C0, 6);
        for (const v of encodeSplat(1, 1, Color(0x000000)).fDc) expect(v).toBeCloseTo(-0.5 / C0, 6);
        const gray = encodeSplat(1, 1, Color(0x808080)).fDc;
        expect(gray[0]).toBeCloseTo((128 / 255 - 0.5) / C0, 6);
    });
});
