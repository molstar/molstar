/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { StructureElement, Link } from 'mol-model/structure';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from 'mol-theme/theme';
import { ScaleLegend } from 'mol-util/color/scale';
import { TableLegend } from 'mol-util/color/tables';
import { getPaletteParams, getPalette } from './util';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives every polymer a unique color based on the position (index) of the polymer in the list of polymers in the structure.'

export const PolymerIndexColorThemeParams = {
    ...getPaletteParams({ scaleList: 'RedYellowBlue' }),
}
export type PolymerIndexColorThemeParams = typeof PolymerIndexColorThemeParams
export function getPolymerIndexColorThemeParams(ctx: ThemeDataContext) {
    return PolymerIndexColorThemeParams // TODO return copy
}
export type PolymerIndexColorThemeProps = PD.Values<typeof PolymerIndexColorThemeParams>

export function PolymerIndexColorTheme(ctx: ThemeDataContext, props: PD.Values<PolymerIndexColorThemeParams>): ColorTheme<PolymerIndexColorThemeParams> {
    let color: LocationColor
    let legend: ScaleLegend | TableLegend | undefined

    if (ctx.structure) {
        const { units } = ctx.structure
        let polymerCount = 0
        for (let i = 0, il = units.length; i <il; ++i) {
            if (units[i].polymerElements.length > 0) ++polymerCount
        }

        const palette = getPalette(polymerCount, props)
        legend = palette.legend
        
        const unitIdColor = new Map<number, Color>()
        for (let i = 0, j = 0, il = units.length; i <il; ++i) {
            if (units[i].polymerElements.length > 0) {
                unitIdColor.set(units[i].id, palette.color(j))
                ++j
            }
        }

        color = (location: Location): Color => {
            let color: Color | undefined
            if (StructureElement.isLocation(location)) {
                color = unitIdColor.get(location.unit.id)
            } else if (Link.isLocation(location)) {
                color = unitIdColor.get(location.aUnit.id)
            }
            return color !== undefined ? color : DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        factory: PolymerIndexColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    }
}

export const PolymerIndexColorThemeProvider: ColorTheme.Provider<PolymerIndexColorThemeParams> = {
    label: 'Polymer Index',
    factory: PolymerIndexColorTheme,
    getParams: getPolymerIndexColorThemeParams,
    defaultValues: PD.getDefaultValues(PolymerIndexColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
}