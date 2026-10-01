/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '@molstar/core/util/color';
import type { Location } from '@molstar/model/model/location';
import { StructureElement, Bond } from '@molstar/model/model/structure';
import type { ColorTheme, LocationColor } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { getPaletteParams, getPalette } from '@molstar/core/util/color/palette';
import type { TableLegend, ScaleLegend } from '@molstar/core/util/legend';
import { ColorLists, getColorListFromName } from '@molstar/core/util/color/lists';
import { ColorThemeCategory } from './categories.js';

const DefaultList = 'dark-2';
const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every chain instance (single chain or collection of single elements) a unique color based on the position (index) of the chain in the list of chains in the structure.';

export const UnitIndexColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type UnitIndexColorThemeParams = typeof UnitIndexColorThemeParams
export function getUnitIndexColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(UnitIndexColorThemeParams);
    if (ctx.structure) {
        if (ctx.structure.root.units.length > ColorLists[DefaultList].list.length) {
            params.palette.defaultValue.name = 'colors';
            params.palette.defaultValue.params = {
                ...params.palette.defaultValue.params,
                list: { kind: 'interpolate', colors: getColorListFromName(DefaultList).list }
            };
        }
    }
    return params;
}

export function UnitIndexColorTheme(ctx: ThemeDataContext, props: PD.Values<UnitIndexColorThemeParams>): ColorTheme<UnitIndexColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    if (ctx.structure) {
        const { units } = ctx.structure.root;
        const palette = getPalette(units.length, props);
        legend = palette.legend;
        const unitIdColor = new Map<number, Color>();
        for (let i = 0, il = units.length; i < il; ++i) {
            unitIdColor.set(units[i].id, palette.color(i));
        }

        color = (location: Location): Color => {
            if (StructureElement.Location.is(location)) {
                return unitIdColor.get(location.unit.id)!;
            } else if (Bond.isLocation(location)) {
                return unitIdColor.get(location.aUnit.id)!;
            }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: UnitIndexColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    };
}

export const UnitIndexColorThemeProvider: ColorTheme.Provider<UnitIndexColorThemeParams, 'unit-index'> = {
    name: 'unit-index',
    label: 'Chain Instance',
    category: ColorThemeCategory.Chain,
    factory: UnitIndexColorTheme,
    getParams: getUnitIndexColorThemeParams,
    defaultValues: PD.getDefaultValues(UnitIndexColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};