/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { StructureElement, Bond, Structure } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { ScaleLegend, TableLegend } from '../../mol-util/legend';
import { Vec3 } from '../../mol-math/linear-algebra';
import { integerDigitCount } from '../../mol-util/number';
import { ColorLists, getColorListFromName } from '../../mol-util/color/lists';

const DefaultList = 'dark-2';
const DefaultColor = Color(0xCCCCCC);
const Description = `Assigns a color based on the operator HKL value of a transformed chain.`;

export const OperatorHklColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type OperatorHklColorThemeParams = typeof OperatorHklColorThemeParams
export function getOperatorHklColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(OperatorHklColorThemeParams);
    if (ctx.structure) {
        if (getOperatorHklSerialMap(ctx.structure.root).map.size > ColorLists[DefaultList].list.length) {
            params.palette.defaultValue.name = 'colors';
            params.palette.defaultValue.params = {
                ...params.palette.defaultValue.params,
                list: { kind: 'interpolate', colors: getColorListFromName(DefaultList).list }
            };
        }
    }
    return params;
}

const hklOffset = 10000;

function hklKey(hkl: Vec3) {
    return hkl.map(v => `${v + hklOffset}`.padStart(5, '0')).join('');
}

function hklKeySplit(key: string) {
    const len = integerDigitCount(hklOffset, 0);
    const h = parseInt(key.substr(0, len));
    const k = parseInt(key.substr(len, len));
    const l = parseInt(key.substr(len + len, len));
    return [ h - hklOffset, k - hklOffset, l - hklOffset ] as Vec3;
}

function formatHkl(hkl: Vec3) {
    return hkl.map(v => v + 5).join('');
}

function getOperatorHklSerialMap(structure: Structure) {
    const map = new Map<string, number>();
    const set = new Set<string>();
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const k = hklKey(structure.units[i].conformation.operator.hkl);
        set.add(k);
    }
    const arr = Array.from(set.values()).sort();
    arr.forEach(k => map.set(k, map.size));
    const min = hklKeySplit(arr[0]);
    const max = hklKeySplit(arr[arr.length - 1]);
    return { min, max, map };
}

export function OperatorHklColorTheme(ctx: ThemeDataContext, props: PD.Values<OperatorHklColorThemeParams>): ColorTheme<OperatorHklColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    if (ctx.structure) {
        const { min, max, map } = getOperatorHklSerialMap(ctx.structure.root);

        const labelTable: string[] = [];
        map.forEach((v, k) => {
            const i = v % map.size;
            const label = formatHkl(hklKeySplit(k));
            if (labelTable[i] === undefined) labelTable[i] = label;
            else labelTable[i] += `, ${label}`;
        });

        props.palette.params.minLabel = formatHkl(min);
        props.palette.params.maxLabel = formatHkl(max);
        props.palette.params.valueLabel = (i: number) => labelTable[i];

        const palette = getPalette(map.size, props);
        legend = palette.legend;

        color = (location: Location): Color => {
            let serial: number | undefined = undefined;
            if (StructureElement.Location.is(location)) {
                const k = hklKey(location.unit.conformation.operator.hkl);
                serial = map.get(k);
            } else if (Bond.isLocation(location)) {
                const k = hklKey(location.aUnit.conformation.operator.hkl);
                serial = map.get(k);
            }
            return serial === undefined ? DefaultColor : palette.color(serial);
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: OperatorHklColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    };
}

export const OperatorHklColorThemeProvider: ColorTheme.Provider<OperatorHklColorThemeParams, 'operator-hkl'> = {
    name: 'operator-hkl',
    label: 'Operator HKL',
    category: ColorTheme.Category.Symmetry,
    factory: OperatorHklColorTheme,
    getParams: getOperatorHklColorThemeParams,
    defaultValues: PD.getDefaultValues(OperatorHklColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};