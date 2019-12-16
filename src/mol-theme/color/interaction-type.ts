/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Bond } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { Color, ColorMap } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { InteractionsProvider } from '../../mol-model-props/computed/interactions';
import { ThemeDataContext } from '../theme';
import { ColorTheme, LocationColor } from '../color';
import { InteractionType } from '../../mol-model-props/computed/interactions/interactions';
import { TableLegend } from '../../mol-util/legend';
import { Task } from '../../mol-task';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Assigns colors according the interaction type of a link.'

const InteractionTypeColors = ColorMap({
    HydrogenBond: 0x2B83BA,
    Hydrophobic: 0x808080,
    HalogenBond: 0x40FFBF,
    Ionic: 0xF0C814,
    MetalCoordination: 0x8C4099,
    CationPi: 0xFF8000,
    PiStacking: 0x8CB366,
    WeakHydrogenBond: 0xC5DDEC,
})

const InteractionTypeColorTable: [string, Color][] = [
    ['Hydrogen Bond', InteractionTypeColors.HydrogenBond],
    ['Hydrophobic', InteractionTypeColors.Hydrophobic],
    ['Halogen Bond', InteractionTypeColors.HalogenBond],
    ['Ionic', InteractionTypeColors.Ionic],
    ['Metal Coordination', InteractionTypeColors.MetalCoordination],
    ['Cation Pi', InteractionTypeColors.CationPi],
    ['Pi Stacking', InteractionTypeColors.PiStacking],
    ['Weak HydrogenBond', InteractionTypeColors.WeakHydrogenBond],
]

function typeColor(type: InteractionType): Color {
    switch (type) {
        case InteractionType.HydrogenBond:
        case InteractionType.WaterHydrogenBond:
        case InteractionType.BackboneHydrogenBond:
            return InteractionTypeColors.HydrogenBond
        case InteractionType.Hydrophobic:
            return InteractionTypeColors.Hydrophobic
        case InteractionType.HalogenBond:
            return InteractionTypeColors.HalogenBond
        case InteractionType.IonicInteraction:
            return InteractionTypeColors.Ionic
        case InteractionType.MetalCoordination:
            return InteractionTypeColors.MetalCoordination
        case InteractionType.CationPi:
            return InteractionTypeColors.CationPi
        case InteractionType.PiStacking:
            return InteractionTypeColors.PiStacking
        case InteractionType.WeakHydrogenBond:
            return InteractionTypeColors.WeakHydrogenBond
        case InteractionType.Unknown:
            return DefaultColor
    }
}

export const InteractionTypeColorThemeParams = { }
export type InteractionTypeColorThemeParams = typeof InteractionTypeColorThemeParams
export function getInteractionTypeColorThemeParams(ctx: ThemeDataContext) {
    return InteractionTypeColorThemeParams // TODO return copy
}

export function InteractionTypeColorTheme(ctx: ThemeDataContext, props: PD.Values<InteractionTypeColorThemeParams>): ColorTheme<InteractionTypeColorThemeParams> {
    let color: LocationColor

    const interactions = ctx.structure ? InteractionsProvider.getValue(ctx.structure) : undefined
    const contextHash = interactions?.version

    if (interactions && interactions.value) {
        const map = interactions.value
        color = (location: Location) => {
            if (Bond.isLocation(location)) {
                const unitInteractions = map.get(location.aUnit.id)
                if (unitInteractions) {
                    const { links, getLinkIndex } = unitInteractions
                    if (links.edgeCount > 0) {
                        const idx = getLinkIndex(location.aIndex, location.bIndex)
                        if (idx !== -1) return typeColor(links.edgeProps.types[idx])
                    }
                }
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        factory: InteractionTypeColorTheme,
        granularity: 'group',
        color: color,
        props: props,
        contextHash,
        description: Description,
        legend: TableLegend(InteractionTypeColorTable)
    }
}

export const InteractionTypeColorThemeProvider: ColorTheme.Provider<InteractionTypeColorThemeParams> = {
    label: 'Interaction Type',
    factory: InteractionTypeColorTheme,
    getParams: getInteractionTypeColorThemeParams,
    defaultValues: PD.getDefaultValues(InteractionTypeColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
    ensureDependencies: (ctx: ThemeDataContext) => {
        return ctx.structure ? InteractionsProvider.attach(ctx.structure.root) : Task.empty()
    }
}