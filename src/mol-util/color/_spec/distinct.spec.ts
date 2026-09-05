/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../color';
import { distinctColors, getDistinctBaseColors, getDistinctGroupColors } from '../distinct';

describe('distinctColors', () => {
    it('returns no colors for non-positive counts', () => {
        expect(distinctColors(0)).toEqual([]);
        expect(distinctColors(-1)).toEqual([]);
    });

    it('generates the requested number of unique colors', () => {
        const colors = distinctColors(4, { clusteringStepCount: 0 });

        expect(colors).toHaveLength(4);
        expect(new Set(colors).size).toBe(4);
    });
});

describe('getDistinctBaseColors', () => {
    const baseColors = [
        0x377eb8, 0xe41a1c, 0x4daf4a, 0x984ea3,
        0xff7f00, 0xffff33, 0xa65628, 0xf781bf,
    ].map(Color);

    it('uses the built-in base palette', () => {
        expect(getDistinctBaseColors(baseColors.length, 0)).toEqual(baseColors);
    });

    it('rotates the base palette by the requested shift', () => {
        expect(getDistinctBaseColors(4, 50)).toEqual([
            baseColors[2], baseColors[3], baseColors[0], baseColors[1],
        ]);
    });

    it('generates additional colors when the built-in palette is too small', () => {
        const colors = getDistinctBaseColors(baseColors.length + 1, 0);

        expect(colors).toHaveLength(baseColors.length + 1);
        expect(new Set(colors).size).toBe(baseColors.length + 1);
    });
});

describe('getDistinctGroupColors', () => {
    it('rotates generated colors by the requested shift', () => {
        const color = Color(0x377eb8);
        const colors = getDistinctGroupColors(4, color, 20, 0);

        expect(getDistinctGroupColors(4, color, 20, 50)).toEqual([
            colors[2], colors[3], colors[0], colors[1],
        ]);
    });

    it('generates grayscale colors for an achromatic base color', () => {
        const colors = getDistinctGroupColors(4, Color(0x808080), 20, 0);

        expect(colors).toHaveLength(4);
        for (const color of colors) {
            const [red, green, blue] = Color.toRgb(color);
            expect(red).toBe(green);
            expect(green).toBe(blue);
        }
    });
});