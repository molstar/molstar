/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from '../theme';
import { ColorListOptions, ColorListName, ColorScale } from 'mol-util/color/scale';
import { StructureElement, Unit } from 'mol-model/structure';
import { missingAccessibleSurfaceAreaValue } from 'mol-model/structure/structure/unit/accessible-surface-area/compute';

const DefaultColor = Color(0xFFFFFF)
const Description = 'Assigns a color based on the relative accessible surface area of a residue.'

export const AccessibleSurfaceAreaColorThemeParams = {
    list: PD.ColorScale<ColorListName>('Rainbow', ColorListOptions),
}
export type AccessibleSurfaceAreaColorThemeParams = typeof AccessibleSurfaceAreaColorThemeParams
export function getAccessibleSurfaceAreaColorThemeParams(ctx: ThemeDataContext) {
    return AccessibleSurfaceAreaColorThemeParams // TODO return copy
}

export function AccessibleSurfaceAreaColorTheme(ctx: ThemeDataContext, props: PD.Values<AccessibleSurfaceAreaColorThemeParams>): ColorTheme<AccessibleSurfaceAreaColorThemeParams> {
    let color: LocationColor
    let scale: ColorScale | undefined = undefined

    if (ctx.structure) {
        scale = ColorScale.create({
            listOrName: props.list,
            minLabel: 'Start',
            maxLabel: 'End'
        })
        const scaleColor = scale.color

        color = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                if (Unit.isAtomic(location.unit)) {
                    const value = location.unit.accessibleSurfaceArea.relativeAccessibleSurfaceArea[location.unit.residueIndex[location.element]];
                    return value !== missingAccessibleSurfaceAreaValue ? scaleColor(value) : DefaultColor;
                }
            }

            return DefaultColor
        }
    } else {
        color = () => DefaultColor
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