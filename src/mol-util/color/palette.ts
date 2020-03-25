/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../param-definition'
import { DistinctColorsParams, distinctColors } from './distinct';
import { ColorScale } from './scale';
import { Color } from '.';
import { ColorListName, ColorListOptionsScale, ColorListOptionsSet, getColorListFromName } from './lists';
import { TableLegend, ScaleLegend } from '../legend';
import { ColorNames } from './names';

type PaletteType = 'generate' | 'scale' | 'set'

const DefaultGetPaletteProps = {
    type: 'generate' as PaletteType,
    scaleList: 'red-yellow-blue' as ColorListName,
    setList: 'set-1' as ColorListName
}
type GetPaletteProps = typeof DefaultGetPaletteProps

const LabelParams = {
    valueLabel: PD.Value((i: number) => `${i + 1}`, { isHidden: true }),
    minLabel: PD.Value('Start', { isHidden: true }),
    maxLabel: PD.Value('End', { isHidden: true }),
}

export function getPaletteParams(props: Partial<GetPaletteProps> = {}) {
    const p = { ...DefaultGetPaletteProps, ...props }
    return {
        palette: PD.MappedStatic(p.type, {
            scale: PD.Group({
                ...LabelParams,
                list: PD.MappedStatic('predefined', {
                    predefined: PD.ColorList<ColorListName>(p.scaleList, ColorListOptionsScale),
                    custom: PD.ObjectList({ color: PD.Color(0x0 as Color) }, ({ color }) => Color.toHexString(color), {
                        defaultValue: [ { color: ColorNames.red }, { color: ColorNames.green }, { color: ColorNames.blue } ]
                    })
                })
            }, { isFlat: true }),
            set: PD.Group({
                ...LabelParams,
                list: PD.ColorList<ColorListName>(p.setList, ColorListOptionsSet),
            }, { isFlat: true }),
            generate: PD.Group({
                ...LabelParams,
                ...DistinctColorsParams,
                maxCount: PD.Numeric(75, { min: 1, max: 250, step: 1 }),
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
        const { list, minLabel, maxLabel } = props.palette.params
        const domain: [number, number] = [0, count - 1]

        let listOrName: ColorListName | Color[];
        if (list.name === 'predefined') {
            listOrName = list.params;
        } else {
            listOrName = list.params.map(c => c.color);
            if (listOrName.length === 0) listOrName = [ColorNames.black];
        }

        const scale = ColorScale.create({ listOrName, domain, minLabel, maxLabel })
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
        const valueLabel = props.palette.params.valueLabel || (i => '' + i);
        const colorsLength = colors.length
        const table: [string, Color][] = []
        for (let i = 0; i < count; ++i) {
            const j = i % colorsLength
            if (table[j] === undefined) {
                table[j] = [valueLabel(i), colors[j]]
            } else {
                table[j][0] += `, ${valueLabel(i)}`
            }
        }
        legend = TableLegend(table)
        color = (i: number) => colors[i % colorsLength]
    }

    return { color, legend }
}