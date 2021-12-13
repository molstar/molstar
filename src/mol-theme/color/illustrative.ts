/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementSymbol } from '../../mol-model/structure/model/types';
import { Color } from '../../mol-util/color';
import { StructureElement, Unit, Bond } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { ColorTheme } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { ChainIdColorTheme, ChainIdColorThemeParams } from './chain-id';

const DefaultIllustrativeColor = Color(0xEEEEEE);
const Description = `Assigns an illustrative color that gives every chain a unique color with lighter carbons (inspired by David Goodsell's Molecule of the Month style).`;

export const IllustrativeColorThemeParams = {
    ...ChainIdColorThemeParams,
    carbonLightness: PD.Numeric(0.8, { min: -6, max: 6, step: 0.1 })
};
export type IllustrativeColorThemeParams = typeof IllustrativeColorThemeParams
export function getIllustrativeColorThemeParams(ctx: ThemeDataContext) {
    return IllustrativeColorThemeParams; // TODO return copy
}

export function IllustrativeColorTheme(ctx: ThemeDataContext, props: PD.Values<IllustrativeColorThemeParams>): ColorTheme<IllustrativeColorThemeParams> {
    const { color: chainIdColor, legend } = ChainIdColorTheme(ctx, props);

    function illustrativeColor(location: Location, typeSymbol: ElementSymbol) {
        const baseColor = chainIdColor(location, false);
        return typeSymbol === 'C' ? Color.lighten(baseColor, props.carbonLightness) : baseColor;
    }

    function color(location: Location): Color {
        if (StructureElement.Location.is(location) && Unit.isAtomic(location.unit)) {
            const typeSymbol = location.unit.model.atomicHierarchy.atoms.type_symbol.value(location.element);
            return illustrativeColor(location, typeSymbol);
        } else if (Bond.isLocation(location) && Unit.isAtomic(location.aUnit)) {
            const elementIndex = location.aUnit.elements[location.aIndex];
            const typeSymbol = location.aUnit.model.atomicHierarchy.atoms.type_symbol.value(elementIndex);
            return illustrativeColor(location, typeSymbol);
        }
        return DefaultIllustrativeColor;
    }

    return {
        factory: IllustrativeColorTheme,
        granularity: 'group',
        preferSmoothing: true,
        color,
        props,
        description: Description,
        legend
    };
}

export const IllustrativeColorThemeProvider: ColorTheme.Provider<IllustrativeColorThemeParams, 'illustrative'> = {
    name: 'illustrative',
    label: 'Illustrative',
    category: ColorTheme.Category.Misc,
    factory: IllustrativeColorTheme,
    getParams: getIllustrativeColorThemeParams,
    defaultValues: PD.getDefaultValues(IllustrativeColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};