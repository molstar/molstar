/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { DistinctColorsParams, distinctColors } from '../../mol-util/color/distinct';
import { ScaleLegend, ColorScale } from '../../mol-util/color/scale';
import { Color } from '../../mol-util/color';
import { TableLegend, ColorListName, ColorListOptionsScale, ColorListOptionsSet, getColorListFromName } from '../../mol-util/color/lists';

const DefaultGetPaletteProps = {
    scaleList: 'red-yellow-blue' as ColorListName,
    setList: 'set-1' as ColorListName
}
type GetPaletteProps = typeof DefaultGetPaletteProps

export function getPaletteParams(props: Partial<GetPaletteProps> = {}) {
    const p = { ...DefaultGetPaletteProps, props }
    return {
        palette: PD.MappedStatic('generate', {
            scale: PD.Group({
                list: PD.ColorList<ColorListName>(p.scaleList, ColorListOptionsScale),
            }, { isFlat: true }),
            set: PD.Group({
                list: PD.ColorList<ColorListName>(p.setList, ColorListOptionsSet),
            }, { isFlat: true }),
            generate: PD.Group({
                ...DistinctColorsParams,
                maxCount: PD.Numeric(75, { min: 1, max: 250, step: 1 })
            }, { isFlat: true })
        }, {
            options: [
                ['scale', 'Interpolate'],
                ['set', 'From Set'],
                ['generate', 'Generate Distinct']
            ]
        })
    }
}

const DefaultPaletteProps = PD.getDefaultValues(getPaletteParams())
type PaletteProps = typeof DefaultPaletteProps

export interface Palette {
    color: (i: number) => Color
    legend?: TableLegend | ScaleLegend
}

export function getPalette(count: number, props: PaletteProps) {
    let color: (i: number) => Color
    let legend: ScaleLegend | TableLegend | undefined

    if (props.palette.name === 'scale') {
        const listOrName = props.palette.params.list
        const domain: [number, number] = [0, count - 1]
        const scale = ColorScale.create({ listOrName, domain, minLabel: 'Start', maxLabel: 'End' })
        legend = scale.legend
        color = scale.color
    } else {
        let colors: Color[]
        if (props.palette.name === 'set') {
            const listOrName = props.palette.params.list
            colors = typeof listOrName === 'string' ? getColorListFromName(listOrName).list : listOrName
        } else {
            count = Math.min(count, props.palette.params.maxCount)
            colors = distinctColors(count, props.palette.params)
        }
        const colorsLength = colors.length
        color = (i: number) => colors[i % colorsLength]
    }

    return { color, legend }
}