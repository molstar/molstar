/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import type { Location } from '@molstar/model/model/location';
import { Color, ColorMap } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { InteractionsProvider } from '@molstar/model/props/computed/interactions';
import type { ThemeDataContext } from '@molstar/graphics/theme/theme';
import type { ColorTheme, LocationColor } from '@molstar/graphics/theme/color';
import { InteractionType } from '@molstar/model/props/computed/interactions/common';
import { TableLegend } from '@molstar/core/util/legend';
import { Interactions, Bridges } from '@molstar/model/props/computed/interactions/interactions';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import { hash2 } from '@molstar/core/data/util';
import { ColorThemeCategory } from '@molstar/graphics/theme/color/categories';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Assigns colors according the interaction type of a link.';

const InteractionTypeColors = ColorMap({
    HydrogenBond: 0x2B83BA,
    Hydrophobic: 0x808080,
    HalogenBond: 0x40FFBF,
    Ionic: 0xF0C814,
    MetalCoordination: 0x8C4099,
    CationPi: 0xFF8000,
    PiStacking: 0x8CB366,
    WeakHydrogenBond: 0xC5DDEC,
    WaterBridge: 0x00CCEE,
});

const InteractionTypeColorTable: [string, Color][] = [
    ['Hydrogen Bond', InteractionTypeColors.HydrogenBond],
    ['Hydrophobic', InteractionTypeColors.Hydrophobic],
    ['Halogen Bond', InteractionTypeColors.HalogenBond],
    ['Ionic', InteractionTypeColors.Ionic],
    ['Metal Coordination', InteractionTypeColors.MetalCoordination],
    ['Cation Pi', InteractionTypeColors.CationPi],
    ['Pi Stacking', InteractionTypeColors.PiStacking],
    ['Weak HydrogenBond', InteractionTypeColors.WeakHydrogenBond],
    ['Water Bridge', InteractionTypeColors.WaterBridge],
];

function typeColor(type: InteractionType): Color {
    switch (type) {
        case InteractionType.HydrogenBond:
            return InteractionTypeColors.HydrogenBond;
        case InteractionType.Hydrophobic:
            return InteractionTypeColors.Hydrophobic;
        case InteractionType.HalogenBond:
            return InteractionTypeColors.HalogenBond;
        case InteractionType.Ionic:
            return InteractionTypeColors.Ionic;
        case InteractionType.MetalCoordination:
            return InteractionTypeColors.MetalCoordination;
        case InteractionType.CationPi:
            return InteractionTypeColors.CationPi;
        case InteractionType.PiStacking:
            return InteractionTypeColors.PiStacking;
        case InteractionType.WeakHydrogenBond:
            return InteractionTypeColors.WeakHydrogenBond;
        case InteractionType.WaterBridge:
            return InteractionTypeColors.WaterBridge;
        case InteractionType.Unknown:
            return DefaultColor;
    }
}

export const InteractionTypeColorThemeParams = { };
export type InteractionTypeColorThemeParams = typeof InteractionTypeColorThemeParams
export function getInteractionTypeColorThemeParams(ctx: ThemeDataContext) {
    return InteractionTypeColorThemeParams; // TODO return copy
}

export function InteractionTypeColorTheme(ctx: ThemeDataContext, props: PD.Values<InteractionTypeColorThemeParams>): ColorTheme<InteractionTypeColorThemeParams> {
    let color: LocationColor;

    const interactions = ctx.structure ? InteractionsProvider.get(ctx.structure) : undefined;
    const contextHash = interactions ? hash2(interactions.id, interactions.version) : -1;

    if (interactions && interactions.value) {
        color = (location: Location) => {
            if (Interactions.isLocation(location)) {
                const { unitsContacts, contacts } = location.data.interactions;
                const { unitA, unitB, indexA, indexB } = location.element;
                if (unitA === unitB) {
                    const links = unitsContacts.get(unitA.id);
                    const idx = links.getDirectedEdgeIndex(indexA, indexB);
                    return typeColor(links.edgeProps.type[idx]);
                } else {
                    const idx = contacts.getEdgeIndex(indexA, unitA.id, indexB, unitB.id);
                    return typeColor(contacts.edges[idx].props.type);
                }
            }
            if (Bridges.isLocation(location)) {
                return typeColor(location.data.bridges[location.element.bridgeIndex].props.type);
            }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: InteractionTypeColorTheme,
        granularity: 'group',
        color: color,
        props: props,
        contextHash,
        description: Description,
        legend: TableLegend(InteractionTypeColorTable)
    };
}

export const InteractionTypeColorThemeProvider: ColorTheme.Provider<InteractionTypeColorThemeParams, 'interaction-type'> = {
    name: 'interaction-type',
    label: 'Interaction Type',
    category: ColorThemeCategory.Misc,
    factory: InteractionTypeColorTheme,
    getParams: getInteractionTypeColorThemeParams,
    defaultValues: PD.getDefaultValues(InteractionTypeColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? InteractionsProvider.attach(ctx, data.structure, void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && InteractionsProvider.ref(data.structure, false)
    }
};