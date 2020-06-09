/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ThemeDataContext } from '../../mol-theme/theme';
import { ColorTheme, LocationColor } from '../../mol-theme/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Color } from '../../mol-util/color';
import { TableLegend } from '../../mol-util/legend';
import { StructureElement, Unit, StructureProperties } from '../../mol-model/structure';
import { isHydrophobic } from '../../mol-model-props/computed/membrane-orientation/ANVIL';
import { Location } from '../../mol-model/location';

const DefaultColor = Color(0x909090);

const HydrophobicColor = Color(0x2166ac);
const HydrophilicColor = Color(0xfee08b);

const ColorLegend = TableLegend([
    ['Data unavailable', DefaultColor],
    ['Hydrophobic', HydrophobicColor],
    ['Hydrophilic', HydrophilicColor],
]);

export function HydrophobicityColorTheme(ctx: ThemeDataContext, props: {}): ColorTheme<{}> {
    let color: LocationColor = () => DefaultColor;
    const { label_comp_id } = StructureProperties.residue;

    if (ctx.structure) {
        color = (location: Location): Color => {
            if (StructureElement.Location.is(location) && Unit.isAtomic(location.unit)) {
                if (isHydrophobic(label_comp_id(location))) {
                    return HydrophobicColor;
                } else {
                    return HydrophilicColor;
                }
            }
            return DefaultColor;
        };
    }

    return {
        factory: HydrophobicityColorTheme,
        granularity: 'group',
        color,
        props,
        description: 'Assigns residue colors according to their hydrophobicity.',
        legend: ColorLegend
    };
}

export const HydrophobicityColorThemeProvider: ColorTheme.Provider<{}, 'hydrophobicity'> = {
    name: 'hydrophobicity',
    label: 'Hydrophobicity',
    category: ColorTheme.Category.Residue,
    factory: HydrophobicityColorTheme,
    getParams: () => ({}),
    defaultValues: PD.getDefaultValues({}),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};