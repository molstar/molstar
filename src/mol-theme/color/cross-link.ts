/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Link } from 'mol-model/structure';

import { Color, ColorScale } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { Vec3 } from 'mol-math/linear-algebra';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from 'mol-theme/theme';
import { ColorListName, ColorListOptions } from 'mol-util/color/scale';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Colors cross-links by the deviation of the observed distance versus the modeled distance (e.g. `ihm_cross_link_restraint.distance_threshold`).'

export const CrossLinkColorThemeParams = {
    domain: PD.Interval([-10, 10]),
    list: PD.Select<ColorListName>('RdYlBu', ColorListOptions),
}
export function getCrossLinkColorThemeParams(ctx: ThemeDataContext) {
    return CrossLinkColorThemeParams // TODO return copy
}
export type CrossLinkColorThemeProps = PD.Values<typeof CrossLinkColorThemeParams>

const distVecA = Vec3.zero(), distVecB = Vec3.zero()
function linkDistance(link: Link.Location) {
    link.aUnit.conformation.position(link.aIndex, distVecA)
    link.bUnit.conformation.position(link.bIndex, distVecB)
    return Vec3.distance(distVecA, distVecB)
}

export function CrossLinkColorTheme(ctx: ThemeDataContext, props: CrossLinkColorThemeProps): ColorTheme<CrossLinkColorThemeProps> {
    let color: LocationColor
    let scale: ColorScale | undefined = undefined

    if (ctx.structure) {
        const crosslinks = ctx.structure.crossLinkRestraints
        scale = ColorScale.create({
            domain: props.domain,
            listOrName: props.list
        })
        const scaleColor = scale.color

        color = (location: Location): Color => {
            if (Link.isLocation(location)) {
                const pairs = crosslinks.getPairs(location.aIndex, location.aUnit, location.bIndex, location.bUnit)
                if (pairs) {
                    return scaleColor(linkDistance(location) - pairs[0].distanceThreshold)
                }
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}

export const CrossLinkColorThemeProvider: ColorTheme.Provider<typeof CrossLinkColorThemeParams> = {
    label: 'Cross Link',
    factory: CrossLinkColorTheme,
    getParams: getCrossLinkColorThemeParams
}