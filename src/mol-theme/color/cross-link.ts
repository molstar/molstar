/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Link } from 'mol-model/structure';

import { Color, ColorScale } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ColorThemeProps, ColorTheme, LocationColor } from '../color';
import { Vec3 } from 'mol-math/linear-algebra';
import { ColorBrewer } from 'mol-util/color/tables';
import { defaults } from 'mol-util';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Colors cross-links by the deviation of the observed distance versus the modeled distance (e.g. `ihm_cross_link_restraint.distance_threshold`).'

const distVecA = Vec3.zero(), distVecB = Vec3.zero()
function linkDistance(link: Link.Location) {
    link.aUnit.conformation.position(link.aIndex, distVecA)
    link.bUnit.conformation.position(link.bIndex, distVecB)
    return Vec3.distance(distVecA, distVecB)
}

export function CrossLinkColorTheme(props: ColorThemeProps): ColorTheme {
    let color: LocationColor
    let scale: ColorScale | undefined = undefined

    if (props.structure) {
        const crosslinks = props.structure.crossLinkRestraints
        scale = ColorScale.create({
            domain: defaults(props.domain, [ -10, 10 ]),
            list: defaults(props.list, ColorBrewer.RdYlBu)
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
        features: { list: true, domain: true, structure: true },
        granularity: 'group',
        color,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}