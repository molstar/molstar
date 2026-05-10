/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec2 } from '../../linear-algebra';
import { pointInPolygon, pointInPolygon2D } from '../polygon';

describe('pointInPolygon', () => {
    it('basic', () => {
        const polygon = [
            -1, -1,
            1, -1,
            1, 1,
            -1, 1
        ];
        expect(pointInPolygon(Vec2.create(0, 0), polygon, 4)).toBe(true);
        expect(pointInPolygon(Vec2.create(2, 2), polygon, 4)).toBe(false);
    });
});

describe('pointInPolygon2D', () => {
    const square: [number, number][] = [[-1, -1], [1, -1], [1, 1], [-1, 1]];

    it('inside square', () => {
        expect(pointInPolygon2D(0, 0, square)).toBe(true);
    });

    it('outside square', () => {
        expect(pointInPolygon2D(2, 2, square)).toBe(false);
    });

    it('degenerate polygon (< 3 vertices)', () => {
        expect(pointInPolygon2D(0, 0, [[0, 0], [1, 0]])).toBe(false);
    });

    it('triangle', () => {
        const tri: [number, number][] = [[0, 0], [1, 0], [0.5, 1]];
        expect(pointInPolygon2D(0.5, 0.4, tri)).toBe(true);
        expect(pointInPolygon2D(0.5, 1.1, tri)).toBe(false);
    });
});