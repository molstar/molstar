/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Color from './color'
import { ColorBrewer } from './tables'

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
            return Color.interpolate(c1, c2, t - tf)
        }
        return {
            color,
            colorToArray: (value: number, array: Helpers.NumberArray, offset: number) => {
                Color.toArray(color(value), array, offset)
            },
            normalizedColorToArray: (value: number, array: Helpers.NumberArray, offset: number) => {
                Color.toArrayNormalized(color(value), array, offset)
            },
        }
    }
}
