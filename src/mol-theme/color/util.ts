/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { DistinctColorsParams, distinctColors } from '../../mol-util/color/distinct';
import { ColorListName, ColorListOptions, ScaleLegend, ColorScale } from '../../mol-util/color/scale';
import { Color } from '../../mol-util/color';
import { TableLegend } from '../../mol-util/color/tables';

const DefaultGetPaletteProps = {
    scaleList: 'RedYellowBlue' as ColorListName
}
type GetPaletteProps = typeof DefaultGetPaletteProps

export function getPaletteParams(props: Partial<GetPaletteProps> = {}) {
    const p = { ...DefaultGetPaletteProps, props }
    return {
        palette: PD.MappedStatic('generate', {
            scale: PD.Group({
                list: PD.ColorScale<ColorListName>(p.scaleList, ColorListOptions),
            }, { isFlat: true }),
            generate: PD.Group({
                ...DistinctColorsParams,
                maxCount: PD.Numeric(75, { min: 1, max: 250, step: 1 })
            }, { isFlat: true })
        }, {
            options: [
                ['scale', 'From Scale'],
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
        count = Math.min(count, props.palette.params.maxCount)
        const colors = distinctColors(count, props.palette.params)
        color = (i: number) => colors[i % count]
    }

    return { color, legend }
}