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
import { Column } from 'mol-data/db';
import { Entities } from 'mol-model/structure/model/properties/common';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives every polymer chain a color based on its `asym_id` value.'

export const PolymerIdColorThemeParams = {
    list: PD.ColorScale<ColorListName>('RedYellowBlue', ColorListOptions),
}
export type PolymerIdColorThemeParams = typeof PolymerIdColorThemeParams
export function getPolymerIdColorThemeParams(ctx: ThemeDataContext) {
    return PolymerIdColorThemeParams // TODO return copy
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

function addPolymerAsymIds(map: Map<string, number>, asymId: Column<string>, entityId: Column<string>, entities: Entities) {
    let j = map.size
    for (let o = 0, ol = asymId.rowCount; o < ol; ++o) {
        const e = entityId.value(o)
        const eI = entities.getEntityIndex(e)
        if (entities.data.type.value(eI) === 'polymer') {
            const k = asymId.value(o)
            if (!map.has(k)) {
                map.set(k, j)
                j += 1
            }
        }
    }
}

export function PolymerIdColorTheme(ctx: ThemeDataContext, props: PD.Values<PolymerIdColorThemeParams>): ColorTheme<PolymerIdColorThemeParams> {
    let color: LocationColor
    const scale = ColorScale.create({ listOrName: props.list, minLabel: 'Start', maxLabel: 'End' })

    if (ctx.structure) {
        // TODO same asym ids in different models should get different color
        const l = StructureElement.create()
        const { models } = ctx.structure
        const polymerAsymIdSerialMap = new Map<string, number>()
        for (let i = 0, il = models.length; i <il; ++i) {
            const m = models[i]
            addPolymerAsymIds(polymerAsymIdSerialMap, m.atomicHierarchy.chains.label_asym_id, m.atomicHierarchy.chains.label_entity_id, m.entities)
            if (m.coarseHierarchy.isDefined) {
                addPolymerAsymIds(polymerAsymIdSerialMap, m.coarseHierarchy.spheres.asym_id, m.coarseHierarchy.spheres.entity_id, m.entities)
                addPolymerAsymIds(polymerAsymIdSerialMap, m.coarseHierarchy.gaussians.asym_id, m.coarseHierarchy.spheres.entity_id, m.entities)
            }
        }
        scale.setDomain(0, polymerAsymIdSerialMap.size - 1)
        const scaleColor = scale.color

        color = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                const asym_id = getAsymId(location.unit)
                return scaleColor(polymerAsymIdSerialMap.get(asym_id(location)) || 0)
            } else if (Link.isLocation(location)) {
                const asym_id = getAsymId(location.aUnit)
                l.unit = location.aUnit
                l.element = location.aUnit.elements[location.aIndex]
                return scaleColor(polymerAsymIdSerialMap.get(asym_id(l)) || 0)
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        factory: PolymerIdColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}

export const PolymerIdColorThemeProvider: ColorTheme.Provider<PolymerIdColorThemeParams> = {
    label: 'Polymer Id',
    factory: PolymerIdColorTheme,
    getParams: getPolymerIdColorThemeParams,
    defaultValues: PD.getDefaultValues(PolymerIdColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
}