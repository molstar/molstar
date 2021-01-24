/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '.';
import { ScaleLegend, TableLegend } from '../legend';
import { ParamDefinition as PD } from '../param-definition';
import { distinctColors, DistinctColorsParams } from './distinct';
import { ColorListName, getColorListFromName } from './lists';
import { ColorScale } from './scale';

type PaletteType = 'generate' | 'colors'

const DefaultGetPaletteProps = {
    type: 'generate' as PaletteType,
    colorList: 'red-yellow-blue' as ColorListName
};
type GetPaletteProps = typeof DefaultGetPaletteProps

export function getPaletteParams(props: Partial<GetPaletteProps> = {}) {
    const p = { ...DefaultGetPaletteProps, ...props };
    return {
        palette: PD.MappedStatic(p.type, {
            colors: PD.Group({
                list: PD.ColorList(p.colorList),
            }, { isFlat: true }),
            generate: PD.Group({
                ...DistinctColorsParams,
                maxCount: PD.Numeric(75, { min: 1, max: 250, step: 1 }),
            }, { isFlat: true })
        }, {
            options: [
                ['colors', 'Color List'],
                ['generate', 'Generate Distinct']
            ]
        })
    };
}

const DefaultPaletteProps = PD.getDefaultValues(getPaletteParams());
type PaletteProps = typeof DefaultPaletteProps

const DefaultLabelOptions = {
    valueLabel: (i: number) => `${i + 1}`,
    minLabel: 'Start',
    maxLabel: 'End'
};
type LabelOptions = typeof DefaultLabelOptions

export interface Palette {
    color: (i: number) => Color
    legend?: TableLegend | ScaleLegend
}

export function getPalette(count: number, props: PaletteProps, labelOptions: Partial<LabelOptions> = {}) {
    let color: (i: number) => Color;
    let legend: ScaleLegend | TableLegend | undefined;

    if (props.palette.name === 'colors' && props.palette.params.list.kind === 'interpolate') {
        const { list } = props.palette.params;
        const domain: [number, number] = [0, count - 1];
        const { minLabel, maxLabel } = { ...DefaultLabelOptions, ...labelOptions };

        let colors = list.colors;
        if (colors.length === 0) colors = getColorListFromName(DefaultGetPaletteProps.colorList).list;

        const scale = ColorScale.create({ listOrName: colors, domain, minLabel, maxLabel });
        legend = scale.legend;
        color = scale.color;
    } else {
        let colors: Color[];
        if (props.palette.name === 'colors') {
            colors = props.palette.params.list.colors.map(c => Array.isArray(c) ? c[0] : c);
            if (colors.length === 0) colors = getColorListFromName('dark-2').list.map(c => Array.isArray(c) ? c[0] : c);
        } else {
            count = Math.min(count, props.palette.params.maxCount);
            colors = distinctColors(count, props.palette.params);
        }
        const valueLabel = labelOptions.valueLabel ?? DefaultLabelOptions.valueLabel;
        const colorsLength = colors.length;
        const table: [string, Color][] = [];
        for (let i = 0; i < count; ++i) {
            const j = i % colorsLength;
            if (table[j] === undefined) {
                table[j] = [valueLabel(i), colors[j]];
            } else {
                table[j][0] += `, ${valueLabel(i)}`;
            }
        }
        legend = TableLegend(table);
        color = (i: number) => colors[i % colorsLength];
    }

    return { color, legend };
}