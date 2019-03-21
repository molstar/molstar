/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ColorTheme } from '../color';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from '../theme';
import { ColorListOptions, ColorListName, ColorScale } from 'mol-util/color/scale';
import { StructureElement, Link, ElementIndex, Unit } from 'mol-model/structure';

const DefaultAccessibleSurfaceAreaColor = Color(0xCCCCCC)
const Description = 'Assigns a color based on the relative accessible surface area of a residue.'

export const AccessibleSurfaceAreaColorThemeParams = {
    list: PD.ColorScale<ColorListName>('Rainbow', ColorListOptions),
}
export type AccessibleSurfaceAreaColorThemeParams = typeof AccessibleSurfaceAreaColorThemeParams
export function getAccessibleSurfaceAreaColorThemeParams(ctx: ThemeDataContext) {
    return AccessibleSurfaceAreaColorThemeParams // TODO return copy
}

export function AccessibleSurfaceAreaColorTheme(ctx: ThemeDataContext, props: PD.Values<AccessibleSurfaceAreaColorThemeParams>): ColorTheme<AccessibleSurfaceAreaColorThemeParams> {
    const scale = ColorScale.create({
        listOrName: props.list,
        minLabel: 'Start',
        maxLabel: 'End',
    })
    const scaleColor = scale.color

    function asaColor(unit: Unit, element: ElementIndex): Color {
        if (Unit.isAtomic(unit)) {
            return scaleColor(unit.model.properties.asa[unit.residueIndex[element]]);
        }
        return DefaultAccessibleSurfaceAreaColor;
    }

    const color = (location: Location): Color => {
        if (StructureElement.isLocation(location)) {
            return asaColor(location.unit, location.element);
        } else if (Link.isLocation(location)) {
            return asaColor(location.aUnit, location.aUnit.elements[location.aIndex]);
        }
        return DefaultAccessibleSurfaceAreaColor;
    }

    return {
        factory: AccessibleSurfaceAreaColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}

export const AccessibleSurfaceAreaColorThemeProvider: ColorTheme.Provider<AccessibleSurfaceAreaColorThemeParams> = {
    label: 'Accessible Surface Area',
    factory: AccessibleSurfaceAreaColorTheme,
    getParams: getAccessibleSurfaceAreaColorThemeParams,
    defaultValues: PD.getDefaultValues(AccessibleSurfaceAreaColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
}