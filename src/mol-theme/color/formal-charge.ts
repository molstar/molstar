/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorScale } from '../../mol-util/color';
import { StructureElement, Unit, Bond, ElementIndex } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import type { ColorTheme } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { ColorThemeCategory } from './categories';
import { ColorLists } from '../../mol-util/color/lists';

const DefaultFormalChargeColor = Color(0xffff99);
const Description = `Assigns a color based on the formal charge of an atom.`;

export const FormalChargeColorThemeParams = {
    domain: PD.Interval([-3, 3]),
    list: PD.ColorList({ kind: 'set', colors: ColorLists['red-white-blue'].list }),
};
export type FormalChargeColorThemeParams = typeof FormalChargeColorThemeParams
export function getFormalChargeColorThemeParams(ctx: ThemeDataContext) {
    return FormalChargeColorThemeParams; // TODO return copy
}

function getFormalCharge(unit: Unit, element: ElementIndex) {
    if (Unit.isAtomic(unit)) {
        return unit.model.atomicHierarchy.atoms.pdbx_formal_charge.value(element);
    } else {
        return 0;
    }
}

export function FormalChargeColorTheme(ctx: ThemeDataContext, props: PD.Values<FormalChargeColorThemeParams>): ColorTheme<FormalChargeColorThemeParams> {
    const scale = ColorScale.create({
        domain: props.domain,
        listOrName: props.list.colors,
    });

    function color(location: Location): Color {
        if (StructureElement.Location.is(location)) {
            const fc = getFormalCharge(location.unit, location.element);
            return fc !== undefined ? scale.color(fc) : DefaultFormalChargeColor;
        } else if (Bond.isLocation(location)) {
            const fc = getFormalCharge(location.aUnit, location.aUnit.elements[location.aIndex]);
            return fc !== undefined ? scale.color(fc) : DefaultFormalChargeColor;
        }
        return DefaultFormalChargeColor;
    }

    return {
        factory: FormalChargeColorTheme,
        granularity: 'group',
        preferSmoothing: true,
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    };
}

export const FormalChargeColorThemeProvider: ColorTheme.Provider<FormalChargeColorThemeParams, 'formal-charge'> = {
    name: 'formal-charge',
    label: 'Formal Charge',
    category: ColorThemeCategory.Atom,
    factory: FormalChargeColorTheme,
    getParams: getFormalChargeColorThemeParams,
    defaultValues: PD.getDefaultValues(FormalChargeColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.models.some(m => m.atomicHierarchy.atoms.pdbx_formal_charge.isDefined)
};