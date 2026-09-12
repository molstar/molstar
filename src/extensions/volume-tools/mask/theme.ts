/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { isPositionLocation } from '../../../mol-geo/util/location-iterator';
import { Location } from '../../../mol-model/location';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ColorThemeCategory } from '../../../mol-theme/color/categories';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { Color } from '../../../mol-util/color';
import { TableLegend } from '../../../mol-util/legend';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { makeLabelAtPosition } from '../voxel-labels';
import { MaskSelection } from './selection';

const Description = 'Colors a volume surface by whether each voxel is selected by the current view polygons.';

export const MaskSelectionColorThemeParams = {
    selectedColor: PD.Color(Color(0xff6b00), { description: 'Color of voxels the polygons select.' }),
    unselectedColor: PD.Color(Color(0xcccccc), { description: 'Color of voxels the polygons leave out; matches the default uniform volume color.' }),
    /** Mirrors `SelectionStore.version`; bumping it forces a color update after the selection changes. */
    version: PD.Numeric(0, {}, { isHidden: true }),
};
export type MaskSelectionColorThemeParams = typeof MaskSelectionColorThemeParams

export function MaskSelectionColorTheme(ctx: ThemeDataContext, props: PD.Values<MaskSelectionColorThemeParams>): ColorTheme<MaskSelectionColorThemeParams> {
    const volume = ctx.volume;
    const store = volume && MaskSelection.get(volume);

    if (!volume || !store) {
        return {
            factory: MaskSelectionColorTheme,
            granularity: 'uniform',
            color: () => props.unselectedColor,
            props,
            description: Description,
        };
    }

    const colors = [props.unselectedColor, props.selectedColor];
    const selectedAt = makeLabelAtPosition(volume, store.selected);

    const color: LocationColor = (location: Location) => {
        if (!isPositionLocation(location)) return props.unselectedColor;
        return colors[selectedAt(location.position)];
    };

    return {
        factory: MaskSelectionColorTheme,
        granularity: 'vertex',
        preferSmoothing: false,
        color,
        props,
        description: Description,
        contextHash: store.version,
        legend: TableLegend([['Selected', props.selectedColor], ['Not selected', props.unselectedColor]]),
    };
}

export const MaskSelectionColorThemeProvider: ColorTheme.Provider<MaskSelectionColorThemeParams, 'mask-selection'> = {
    name: 'mask-selection',
    label: 'Mask Selection',
    category: ColorThemeCategory.Misc,
    factory: MaskSelectionColorTheme,
    getParams: () => MaskSelectionColorThemeParams,
    defaultValues: PD.getDefaultValues(MaskSelectionColorThemeParams),
    // Applicable to any volume so the theme survives param normalisation before a selection exists.
    isApplicable: (ctx: ThemeDataContext) => !!ctx.volume,
};
