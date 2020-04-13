/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

const LabelParams = {
    valueLabel: PD.Value((i: number) => `${i + 1}`, { isHidden: true }),
    minLabel: PD.Value('Start', { isHidden: true }),
    maxLabel: PD.Value('End', { isHidden: true })
};

export function getPaletteParams(props: Partial<GetPaletteProps> = {}) {
    const p = { ...DefaultGetPaletteProps, ...props };
    return {
        palette: PD.MappedStatic(p.type, {
            colors: PD.Group({
                ...LabelParams,
                list: PD.ColorList(p.colorList),
            }, { isFlat: true }),
            generate: PD.Group({
                ...LabelParams,
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

export interface Palette {
    color: (i: number) => Color
    legend?: TableLegend | ScaleLegend
}

export function getPalette(count: number, props: PaletteProps) {
    let color: (i: number) => Color;
    let legend: ScaleLegend | TableLegend | undefined;

    if (props.palette.name === 'colors' && props.palette.params.list.kind === 'interpolate') {
        const { list, minLabel, maxLabel } = props.palette.params;
        const domain: [number, number] = [0, count - 1];

        let colors = list.colors;
        if (colors.length === 0) colors = getColorListFromName(DefaultGetPaletteProps.colorList).list;

        const scale = ColorScale.create({ listOrName: colors, domain, minLabel, maxLabel });
        legend = scale.legend;
        color = scale.color;
    } else {
        let colors: Color[];
        if (props.palette.name === 'colors') {
            colors = props.palette.params.list.colors;
            if (colors.length === 0) colors = getColorListFromName('dark-2').list;
        } else {
            count = Math.min(count, props.palette.params.maxCount);
            colors = distinctColors(count, props.palette.params);
        }
        const valueLabel = props.palette.params.valueLabel || (i => '' + i);
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