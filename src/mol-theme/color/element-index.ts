/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorScale, Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { StructureElement, Link } from '../../mol-model/structure';
import { OrderedSet } from '../../mol-data/int';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { ThemeDataContext } from '../../mol-theme/theme';
import { ColorListOptions, ColorListName } from '../../mol-util/color/lists';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives every element (atom or coarse sphere/gaussian) a unique color based on the position (index) of the element in the list of elements in the structure.'

export const ElementIndexColorThemeParams = {
    list: PD.ColorScale<ColorListName>('red-yellow-blue', ColorListOptions),
}
export type ElementIndexColorThemeParams = typeof ElementIndexColorThemeParams
export function getElementIndexColorThemeParams(ctx: ThemeDataContext) {
    return ElementIndexColorThemeParams // TODO return copy
}

export function ElementIndexColorTheme(ctx: ThemeDataContext, props: PD.Values<ElementIndexColorThemeParams>): ColorTheme<ElementIndexColorThemeParams> {
    let color: LocationColor
    let scale = ColorScale.create({ listOrName: props.list, minLabel: 'Start', maxLabel: 'End' })

    if (ctx.structure) {
        const { units } = ctx.structure
        const unitCount = units.length
        const cummulativeElementCount = new Map<number, number>()
        const unitIdIndex = new Map<number, number>()

        let elementCount = 0
        for (let i = 0; i < unitCount; ++i) {
            cummulativeElementCount.set(i, elementCount)
            elementCount += units[i].elements.length
            unitIdIndex.set(units[i].id, i)
        }
        scale.setDomain(0, elementCount - 1)
        const scaleColor = scale.color

        color = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                const unitIndex = unitIdIndex.get(location.unit.id)!
                const unitElementIndex = OrderedSet.findPredecessorIndex(location.unit.elements, location.element)
                return scaleColor(cummulativeElementCount.get(unitIndex)! + unitElementIndex)
            } else if (Link.isLocation(location)) {
                const unitIndex = unitIdIndex.get(location.aUnit.id)!
                return scaleColor(cummulativeElementCount.get(unitIndex)! + location.aIndex)
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        factory: ElementIndexColorTheme,
        granularity: 'groupInstance',
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}

export const ElementIndexColorThemeProvider: ColorTheme.Provider<ElementIndexColorThemeParams> = {
    label: 'Element Index',
    factory: ElementIndexColorTheme,
    getParams: getElementIndexColorThemeParams,
    defaultValues: PD.getDefaultValues(ElementIndexColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
}