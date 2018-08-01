/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Create points for a polygon:
 * 3 for a triangle, 4 for a rectangle, 5 for a pentagon, 6 for a hexagon...
 */
export function polygon(sideCount: number, shift: boolean) {
    const points = new Float32Array(sideCount * 2)
    const radius = sideCount <= 4 ? Math.sqrt(2) / 2 : 0.6

    const offset = shift ? 0 : 1

    for (let i = 0, il = 2 * sideCount; i < il; i += 2) {
        const c = (i + offset) / sideCount * Math.PI
        points[i] = Math.cos(c) * radius
        points[i + 1] = Math.sin(c) * radius
    }
    return points
}