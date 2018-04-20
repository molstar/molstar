/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** RGB color triplet expressed as a single number */
type Color = number

namespace Color {
    export function toRgb(hexColor: Color) {
        return [ hexColor >> 16 & 255, hexColor >> 8 & 255, hexColor & 255 ]
    }

    export function toRgbNormalized(hexColor: Color) {
        return [ (hexColor >> 16 & 255) / 255, (hexColor >> 8 & 255) / 255, (hexColor & 255) / 255 ]
    }

    export function fromRgb(r: number, g: number, b: number ): Color {
        return (r << 16) | (g << 8) | b
    }

    /** Copies hex color to rgb array */
    export function toArray(hexColor: Color, array: Helpers.NumberArray, offset: number) {
        array[ offset ] = (hexColor >> 16 & 255)
        array[ offset + 1 ] = (hexColor >> 8 & 255)
        array[ offset + 2 ] = (hexColor & 255)
    }

    /** Copies normalized (0 to 1) hex color to rgb array */
    export function toArrayNormalized(hexColor: Color, array: Helpers.NumberArray, offset: number) {
        array[ offset ] = (hexColor >> 16 & 255) / 255
        array[ offset + 1 ] = (hexColor >> 8 & 255) / 255
        array[ offset + 2 ] = (hexColor & 255) / 255
    }

    /** Linear interpolation between two colors */
    export function interpolate(c1: Color, c2: Color, t: number): Color {
        const r1 = c1 >> 16 & 255
        const g1 = c1 >> 8 & 255
        const b1 = c1 & 255
        const r2 = c2 >> 16 & 255
        const g2 = c2 >> 8 & 255
        const b2 = c2 & 255

        const r = r1 + (r2 - r1) * t
        const g = g1 + (g2 - g1) * t
        const b = b1 + (b2 - b1) * t

        return (r << 16) | (g << 8) | b
    }
}

export default Color