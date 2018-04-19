/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorBrewer } from './tables'

export { ElementColor } from './structure/element';

export function colorToRgb(hexColor: number) {
    return { r: hexColor >> 16 & 255, g: hexColor >> 8 & 255, b: hexColor & 255 }
}

/** Copies hex color to rgb array */
export function colorToArray(hexColor: number, array: Helpers.NumberArray, offset: number) {

    array[ offset ] = (hexColor >> 16 & 255)
    array[ offset + 1 ] = (hexColor >> 8 & 255)
    array[ offset + 2 ] = (hexColor & 255)
}

/** Copies normalized (0 to 1) hex color to rgb array */
export function normalizedColorToArray(hexColor: number, array: Helpers.NumberArray, offset: number) {
    array[ offset ] = (hexColor >> 16 & 255) / 255
    array[ offset + 1 ] = (hexColor >> 8 & 255) / 255
    array[ offset + 2 ] = (hexColor & 255) / 255
}

/** Linear interpolation between two colors */
export function interpolate(c1: number, c2: number, t: number) {
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

export type Color = number

export interface ColorScale {
    /** Returns hex color for given value */
    color: (value: number) => Color
    /** Copies color to rgb int8 array */
    colorToArray: (value: number, array: Helpers.NumberArray, offset: number) => void
    /** Copies normalized (0 to 1) hex color to rgb array */
    normalizedColorToArray: (value: number, array: Helpers.NumberArray, offset: number) => void
}

export const DefaultColorScale = {
    domain: [0, 1],
    reverse: false,
    colors: ColorBrewer.RdYlBu as Color[]
}
export type ColorScaleProps = Partial<typeof DefaultColorScale>

export namespace ColorScale {
    export function create(props: ColorScaleProps): ColorScale {
        const { domain, reverse, colors } = { ...DefaultColorScale, ...props }
        const [ min, max ] = reverse ? domain.slice().reverse() : domain
        const count1 = colors.length - 1

        function color(value: number) {
            const t = ((value - min) / (max - min)) * count1
            const tf = Math.floor(t)
            const c1 = colors[tf]
            const c2 = colors[Math.ceil(t)]
            return interpolate(c1, c2, t - tf)
        }
        return {
            color,
            colorToArray: (value: number, array: Helpers.NumberArray, offset: number) => {
                colorToArray(color(value), array, offset)
            },
            normalizedColorToArray: (value: number, array: Helpers.NumberArray, offset: number) => {
                normalizedColorToArray(color(value), array, offset)
            },
        }
    }
}
