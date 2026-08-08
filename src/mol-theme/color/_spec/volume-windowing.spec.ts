/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { applyVolumeWindowing, getVolumeWindowingPalette } from '../volume-windowing';
import { Color } from '../../../mol-util/color';

describe('volume-windowing theme', () => {
    it('maps the full range linearly by default', () => {
        expect(applyVolumeWindowing(0, [0, 255], 1, false)).toBeCloseTo(0, 6);
        expect(applyVolumeWindowing(0.5, [0, 255], 1, false)).toBeCloseTo(0.5, 6);
        expect(applyVolumeWindowing(1, [0, 255], 1, false)).toBeCloseTo(1, 6);
    });

    it('applies black and white clipping in IMOD-style display levels', () => {
        expect(applyVolumeWindowing(32 / 255, [64, 192], 1, false)).toBeCloseTo(0, 6);
        expect(applyVolumeWindowing(64 / 255, [64, 192], 1, false)).toBeCloseTo(0, 6);
        expect(applyVolumeWindowing(128 / 255, [64, 192], 1, false)).toBeCloseTo(0.5, 6);
        expect(applyVolumeWindowing(192 / 255, [64, 192], 1, false)).toBeCloseTo(1, 6);
        expect(applyVolumeWindowing(224 / 255, [64, 192], 1, false)).toBeCloseTo(1, 6);
    });

    it('applies gamma and invert after level mapping', () => {
        expect(applyVolumeWindowing(0.25, [0, 255], 2, false)).toBeCloseTo(0.5, 6);
        expect(applyVolumeWindowing(0.25, [0, 255], 2, true)).toBeCloseTo(0.5, 6);
        expect(applyVolumeWindowing(0, [0, 255], 1, true)).toBeCloseTo(1, 6);
        expect(applyVolumeWindowing(1, [0, 255], 1, true)).toBeCloseTo(0, 6);
    });

    it('creates a grayscale palette that respects invert', () => {
        const normal = getVolumeWindowingPalette([0, 255], 1, false);
        const inverted = getVolumeWindowingPalette([0, 255], 1, true);

        expect(normal).toHaveLength(256);
        expect(normal[0]).toBe(Color(0x000000));
        expect(normal[255]).toBe(Color(0xFFFFFF));
        expect(inverted[0]).toBe(Color(0xFFFFFF));
        expect(inverted[255]).toBe(Color(0x000000));
    });
});
