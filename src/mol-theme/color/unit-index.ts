/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorScale, Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { StructureElement, Link } from 'mol-model/structure';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from 'mol-theme/theme';
import { ColorListOptions, ColorListName } from 'mol-util/color/scale';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives every unit (single chain or collection of single elements) a unique color based on the position (index) of the unit in the list of units in the structure.'

export const UnitIndexColorThemeParams = {
    list: PD.Select<ColorListName>('RdYlBu', ColorListOptions),
}
export type UnitIndexColorThemeParams = typeof UnitIndexColorThemeParams
export function getUnitIndexColorThemeParams(ctx: ThemeDataContext) {
    return UnitIndexColorThemeParams // TODO return copy
}

export function UnitIndexColorTheme(ctx: ThemeDataContext, props: PD.Values<UnitIndexColorThemeParams>): ColorTheme<UnitIndexColorThemeParams> {
    let color: LocationColor
    const scale = ColorScale.create({ listOrName: props.list, minLabel: 'Start', maxLabel: 'End' })

    if (ctx.structure) {
        const { units } = ctx.structure
        scale.setDomain(0, units.length - 1)
        const unitIdColor = new Map<number, Color>()
        for (let i = 0, il = units.length; i <il; ++i) {
            unitIdColor.set(units[i].id, scale.color(i))
        }

        color = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                return unitIdColor.get(location.unit.id)!
            } else if (Link.isLocation(location)) {
                return unitIdColor.get(location.aUnit.id)!
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        factory: UnitIndexColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}

export const UnitIndexColorThemeProvider: ColorTheme.Provider<UnitIndexColorThemeParams> = {
    label: 'Unit Index',
    factory: UnitIndexColorTheme,
    getParams: getUnitIndexColorThemeParams,
    defaultValues: PD.getDefaultValues(UnitIndexColorThemeParams)
}