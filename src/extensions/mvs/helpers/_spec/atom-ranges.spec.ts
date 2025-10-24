/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { ElementRanges } from '../element-ranges';


describe('union', () => {
    it('union non-overlapping', async () => {
        const a = {
            from: [0, 20, 40, 60, 80],
            to: [10, 30, 50, 70, 90],
        } as ElementRanges;
        const b = {
            from: [11, 37, 51, 205],
            to: [15, 39, 55, 210],
        } as ElementRanges;
        const c = {
            from: [-10, 200, 300],
            to: [-5, 202, 305],
        } as ElementRanges;
        const result = {
            from: [-10, 0, 11, 20, 37, 40, 51, 60, 80, 200, 205, 300],
            to: [-5, 10, 15, 30, 39, 50, 55, 70, 90, 202, 210, 305],
        } as ElementRanges;
        expect(ElementRanges.union([a, b, c])).toEqual(result);
    });
    it('union overlapping', async () => {
        const a = {
            from: [0, 20, 40, 60, 80],
            to: [10, 30, 50, 70, 90],
        } as ElementRanges;
        const b = {
            from: [10, 37, 51, 84, 205],
            to: [15, 40, 55, 88, 220],
        } as ElementRanges;
        const c = {
            from: [-10, 67, 200, 300],
            to: [5, 80, 210, 305],
        } as ElementRanges;
        const result = {
            from: [-10, 20, 37, 51, 60, 200, 300],
            to: [15, 30, 50, 55, 90, 220, 305],
        } as ElementRanges;
        expect(ElementRanges.union([a, b, c])).toEqual(result);
    });
});

