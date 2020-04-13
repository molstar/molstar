/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Location } from '../../../mol-model/location';
import { Color, ColorMap } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { InteractionsProvider } from '../interactions';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { InteractionType } from '../interactions/common';
import { TableLegend } from '../../../mol-util/legend';
import { Interactions } from '../interactions/interactions';
import { CustomProperty } from '../../common/custom-property';

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
    const contextHash = interactions?.version;

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
                    const idx = contacts.getEdgeIndex(indexA, unitA, indexB, unitB);
                    return typeColor(contacts.edges[idx].props.type);
                }
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
    category: ColorTheme.Category.Misc,
    factory: InteractionTypeColorTheme,
    getParams: getInteractionTypeColorThemeParams,
    defaultValues: PD.getDefaultValues(InteractionTypeColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? InteractionsProvider.attach(ctx, data.structure, void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && data.structure.customPropertyDescriptors.reference(InteractionsProvider.descriptor, false)
    }
};