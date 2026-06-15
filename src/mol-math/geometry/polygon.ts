/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { NumberArray } from '../../mol-util/type-helpers';
import { Vec2 } from '../linear-algebra';

/** raycast along x-axis and apply even-odd rule */
export function pointInPolygon(point: Vec2, polygon: NumberArray, count: number): boolean {
    const [x, y] = point;
    let inside = false;

    for (let i = 0, j = count - 1; i < count; j = i++) {
        const xi = polygon[i * 2], yi = polygon[i * 2 + 1];
        const xj = polygon[j * 2], yj = polygon[j * 2 + 1];

        if (((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
            inside = !inside;
        }
    }
    return inside;
}

/** raycast along x-axis and apply even-odd rule; accepts an array of [x, y] tuples */
export function pointInPolygon2D(px: number, py: number, polygon: readonly [number, number][]): boolean {
    const n = polygon.length;
    if (n < 3) return false;
    let inside = false;
    for (let i = 0, j = n - 1; i < n; j = i++) {
        const xi = polygon[i][0], yi = polygon[i][1];
        const xj = polygon[j][0], yj = polygon[j][1];
        if (((yi > py) !== (yj > py)) && (px < ((xj - xi) * (py - yi)) / (yj - yi) + xi))
            inside = !inside;
    }
    return inside;
}
