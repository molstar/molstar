/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { UniformColorTheme, UniformColorThemeParams } from './uniform';
import { assertUnreachable } from '../../mol-util/type-helpers';
import { EntityIdColorTheme, EntityIdColorThemeParams } from './entity-id';
import { MoleculeTypeColorTheme, MoleculeTypeColorThemeParams } from './molecule-type';
import { EntitySourceColorTheme, EntitySourceColorThemeParams } from './entity-source';

const DefaultIllustrativeColor = Color(0xEEEEEE);
const Description = `Assigns an illustrative color that gives every chain a color based on the choosen style but with lighter carbons (inspired by David Goodsell's Molecule of the Month style).`;

export const IllustrativeColorThemeParams = {
    style: PD.MappedStatic('entity-id', {
        uniform: PD.Group(UniformColorThemeParams),
        'chain-id': PD.Group(ChainIdColorThemeParams),
        'entity-id': PD.Group(EntityIdColorThemeParams),
        'entity-source': PD.Group(EntitySourceColorThemeParams),
        'molecule-type': PD.Group(MoleculeTypeColorThemeParams),
    }),
    carbonLightness: PD.Numeric(0.8, { min: -6, max: 6, step: 0.1 })
};
export type IllustrativeColorThemeParams = typeof IllustrativeColorThemeParams
export function getIllustrativeColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(IllustrativeColorThemeParams);
    return params;
}

export function IllustrativeColorTheme(ctx: ThemeDataContext, props: PD.Values<IllustrativeColorThemeParams>): ColorTheme<IllustrativeColorThemeParams> {
    const { color: styleColor, legend } =
        props.style.name === 'uniform' ? UniformColorTheme(ctx, props.style.params) :
            props.style.name === 'chain-id' ? ChainIdColorTheme(ctx, props.style.params) :
                props.style.name === 'entity-id' ? EntityIdColorTheme(ctx, props.style.params) :
                    props.style.name === 'entity-source' ? EntitySourceColorTheme(ctx, props.style.params) :
                        props.style.name === 'molecule-type' ? MoleculeTypeColorTheme(ctx, props.style.params) :
                            assertUnreachable(props.style);

    function illustrativeColor(location: Location, typeSymbol: ElementSymbol) {
        const baseColor = styleColor(location, false);
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