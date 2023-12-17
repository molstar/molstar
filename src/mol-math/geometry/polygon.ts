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
