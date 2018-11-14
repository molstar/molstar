/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureProperties, StructureElement, Link } from 'mol-model/structure';

import { ColorScale, Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from 'mol-theme/theme';
import { ColorListOptions, ColorListName } from 'mol-util/color/scale';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives every chain a color based on its `asym_id` value.'

export const ChainIdColorThemeParams = {
    list: PD.Select<ColorListName>('Color Scale', '', 'RdYlBu', ColorListOptions),
}
export function getChainIdColorThemeParams(ctx: ThemeDataContext) {
    return ChainIdColorThemeParams // TODO return copy
}
export type ChainIdColorThemeProps = PD.DefaultValues<typeof ChainIdColorThemeParams>

function getAsymId(unit: Unit): StructureElement.Property<string> {
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            return StructureProperties.chain.label_asym_id
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return StructureProperties.coarse.asym_id
    }
}

export function ChainIdColorTheme(ctx: ThemeDataContext, props: ChainIdColorThemeProps): ColorTheme<ChainIdColorThemeProps> {
    let color: LocationColor
    const scale = ColorScale.create({ listOrName: props.list, minLabel: 'Start', maxLabel: 'End' })

    if (ctx.structure) {
        const l = StructureElement.create()
        const { models } = ctx.structure
        const asymIdSerialMap = new Map<string, number>()
        let j = 0
        for (let i = 0, il = models.length; i <il; ++i) {
            models[i].properties.asymIdSerialMap.forEach((v, k) => {
                if (!asymIdSerialMap.has(k)) {
                    asymIdSerialMap.set(k, j)
                    j += 1
                }
            })
        }
        scale.setDomain(0, asymIdSerialMap.size - 1)
        const scaleColor = scale.color

        color = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                const asym_id = getAsymId(location.unit)
                return scaleColor(asymIdSerialMap.get(asym_id(location)) || 0)
            } else if (Link.isLocation(location)) {
                const asym_id = getAsymId(location.aUnit)
                l.unit = location.aUnit
                l.element = location.aUnit.elements[location.aIndex]
                return scaleColor(asymIdSerialMap.get(asym_id(l)) || 0)
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

export const ChainIdColorThemeProvider: ColorTheme.Provider<typeof ChainIdColorThemeParams> = {
    factory: ChainIdColorTheme, params: getChainIdColorThemeParams
}