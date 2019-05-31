/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureProperties, StructureElement, Link } from 'mol-model/structure';

import { Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from 'mol-theme/theme';
import { ScaleLegend } from 'mol-util/color/scale';
import { Column } from 'mol-data/db';
import { getPaletteParams, getPalette } from './util';
import { TableLegend } from 'mol-util/color/tables';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives every chain a color based on its `asym_id` value.'

export const ChainIdColorThemeParams = {
    ...getPaletteParams({ scaleList: 'RedYellowBlue' }),
}
export type ChainIdColorThemeParams = typeof ChainIdColorThemeParams
export function getChainIdColorThemeParams(ctx: ThemeDataContext) {
    return ChainIdColorThemeParams // TODO return copy
}

function getAsymId(unit: Unit): StructureElement.Property<string> {
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            return StructureProperties.chain.label_asym_id
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return StructureProperties.coarse.asym_id
    }
}

function addAsymIds(map: Map<string, number>, data: Column<string>) {
    let j = map.size
    for (let o = 0, ol = data.rowCount; o < ol; ++o) {
        const k = data.value(o)
        if (!map.has(k)) {
            map.set(k, j)
            j += 1
        }
    }
}

export function ChainIdColorTheme(ctx: ThemeDataContext, props: PD.Values<ChainIdColorThemeParams>): ColorTheme<ChainIdColorThemeParams> {
    let color: LocationColor
    let legend: ScaleLegend | TableLegend | undefined

    if (ctx.structure) {
        // TODO same asym ids in different models should get different color
        const l = StructureElement.create()
        const { models } = ctx.structure
        const asymIdSerialMap = new Map<string, number>()
        for (let i = 0, il = models.length; i <il; ++i) {
            const m = models[i]
            addAsymIds(asymIdSerialMap, m.atomicHierarchy.chains.label_asym_id)
            if (m.coarseHierarchy.isDefined) {
                addAsymIds(asymIdSerialMap, m.coarseHierarchy.spheres.asym_id)
                addAsymIds(asymIdSerialMap, m.coarseHierarchy.gaussians.asym_id)
            }
        }

        const palette = getPalette(asymIdSerialMap.size, props)
        legend = palette.legend

        color = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                const asym_id = getAsymId(location.unit)
                return palette.color(asymIdSerialMap.get(asym_id(location)) || 0)
            } else if (Link.isLocation(location)) {
                const asym_id = getAsymId(location.aUnit)
                l.unit = location.aUnit
                l.element = location.aUnit.elements[location.aIndex]
                return palette.color(asymIdSerialMap.get(asym_id(l)) || 0)
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        factory: ChainIdColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend
    }
}

export const ChainIdColorThemeProvider: ColorTheme.Provider<ChainIdColorThemeParams> = {
    label: 'Chain Id',
    factory: ChainIdColorTheme,
    getParams: getChainIdColorThemeParams,
    defaultValues: PD.getDefaultValues(ChainIdColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
}