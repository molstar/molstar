/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { StructureElement, Bond, Structure } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import type { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { ScaleLegend, TableLegend } from '../../mol-util/legend';
import { ColorThemeCategory } from './categories';

const DefaultList = 'many-distinct';
const DefaultColor = Color(0xCCCCCC);
const Description = `Assigns a color based on the operator name of a transformed chain.`;

export const OperatorNameColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type OperatorNameColorThemeParams = typeof OperatorNameColorThemeParams
export function getOperatorNameColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(OperatorNameColorThemeParams);
    return params;
}

function getOperatorNameSerialMap(structure: Structure) {
    const map = new Map<string, number>();
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const name = structure.units[i].conformation.operator.name;
        if (!map.has(name)) map.set(name, map.size);
    }
    return map;
}

export function OperatorNameColorTheme(ctx: ThemeDataContext, props: PD.Values<OperatorNameColorThemeParams>): ColorTheme<OperatorNameColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    if (ctx.structure) {
        const operatorNameSerialMap = getOperatorNameSerialMap(ctx.structure.root);

        const labelTable = Array.from(operatorNameSerialMap.keys());
        const valueLabel = (i: number) => labelTable[i];

        const palette = getPalette(operatorNameSerialMap.size, props, { valueLabel });
        legend = palette.legend;

        color = (location: Location): Color => {
            let serial: number | undefined = undefined;
            if (StructureElement.Location.is(location)) {
                const name = location.unit.conformation.operator.name;
                serial = operatorNameSerialMap.get(name);
            } else if (Bond.isLocation(location)) {
                const name = location.aUnit.conformation.operator.name;
                serial = operatorNameSerialMap.get(name);
            }
            return serial === undefined ? DefaultColor : palette.color(serial);
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: OperatorNameColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    };
}

export const OperatorNameColorThemeProvider: ColorTheme.Provider<OperatorNameColorThemeParams, 'operator-name'> = {
    name: 'operator-name',
    label: 'Operator Name',
    category: ColorThemeCategory.Symmetry,
    factory: OperatorNameColorTheme,
    getParams: getOperatorNameColorThemeParams,
    defaultValues: PD.getDefaultValues(OperatorNameColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};