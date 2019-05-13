/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from './color'
import { ColorBrewer, ColorMatplotlib, ColorOther } from './tables'
import { defaults } from 'mol-util';
import { NumberArray } from 'mol-util/type-helpers';

export type ColorListName = (
    keyof typeof ColorBrewer | keyof typeof ColorMatplotlib | keyof typeof ColorOther
)
export const ColorListNames = [
    ...Object.keys(ColorBrewer), ...Object.keys(ColorMatplotlib), ...Object.keys(ColorOther)
]
export const ColorListOptions = ColorListNames.map(n => [n, n] as [ColorListName, string])

export function getColorListFromName(name: ColorListName) {
    if (name in ColorBrewer) {
        return ColorBrewer[name as keyof typeof ColorBrewer]
    } else if (name in ColorMatplotlib) {
        return ColorMatplotlib[name as keyof typeof ColorMatplotlib]
    } else if (name in ColorOther) {
        return ColorOther[name as keyof typeof ColorOther]
    }
    console.warn(`unknown color list named '${name}'`)
    return ColorBrewer.RedYellowBlue
}

//

export interface ScaleLegend {
    kind: 'scale-legend'
    minLabel: string,
    maxLabel: string,
    colors: Color[]
}
export function ScaleLegend(minLabel: string, maxLabel: string, colors: Color[]): ScaleLegend {
    return { kind: 'scale-legend', minLabel, maxLabel, colors }
}

export interface ColorScale {
    /** Returns hex color for given value */
    color: (value: number) => Color
    /** Copies color to rgb int8 array */
    colorToArray: (value: number, array: NumberArray, offset: number) => void
    /** Copies normalized (0 to 1) hex color to rgb array */
    normalizedColorToArray: (value: number, array: NumberArray, offset: number) => void
    /**  */
    setDomain: (min: number, max: number) => void
    /** Legend */
    readonly legend: ScaleLegend
}

export const DefaultColorScaleProps = {
    domain: [0, 1] as [number, number],
    reverse: false,
    listOrName: ColorBrewer.RedYellowBlue as Color[] | ColorListName,
    minLabel: '' as string | undefined,
    maxLabel: '' as string | undefined,
}
export type ColorScaleProps = Partial<typeof DefaultColorScaleProps>

export namespace ColorScale {
    export function create(props: ColorScaleProps): ColorScale {
        const { domain, reverse, listOrName } = { ...DefaultColorScaleProps, ...props }
        const list = typeof listOrName === 'string' ? getColorListFromName(listOrName) : listOrName

        const colors = reverse ? list.slice().reverse() : list
        const count1 = colors.length - 1

        let diff = 0, min = 0, max = 0
        function setDomain(_min: number, _max: number) {
            min = _min
            max = _max
            diff = (max - min) || 1
        }
        setDomain(domain[0], domain[1])

        const minLabel = defaults(props.minLabel, min.toString())
        const maxLabel = defaults(props.maxLabel, max.toString())

        function color(value: number) {
            const t = Math.min(colors.length - 1, Math.max(0, ((value - min) / diff) * count1))
            const tf = Math.floor(t)
            const c1 = colors[tf]
            const c2 = colors[Math.ceil(t)]
            return Color.interpolate(c1, c2, t - tf)
        }
        return {
            color,
            colorToArray: (value: number, array: NumberArray, offset: number) => {
                Color.toArray(color(value), array, offset)
            },
            normalizedColorToArray: (value: number, array: NumberArray, offset: number) => {
                Color.toArrayNormalized(color(value), array, offset)
            },
            setDomain,
            get legend() { return ScaleLegend(minLabel, maxLabel, colors) }
        }
    }
}
