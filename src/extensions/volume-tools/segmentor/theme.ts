/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ColorThemeCategory } from '../../../mol-theme/color/categories';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Color } from '../../../mol-util/color';
import { TableLegend } from '../../../mol-util/legend';
import { Location } from '../../../mol-model/location';
import { isPositionLocation } from '../../../mol-geo/util/location-iterator';
import { makeLabelAtPosition } from '../voxel-labels';
import { BodyLabels } from './labels';
import { MaxBodyId } from './types';

const Description = 'Colors a volume surface by the body each voxel is assigned to.';

export const BodyLabelColorThemeParams = {
    unassignedColor: PD.Color(Color(0x9a9a9a), { description: 'Color of voxels not assigned to any body.' }),
    /** Mirrors `LabelStore.version`; bumping it forces a color update after labels change. */
    version: PD.Numeric(0, {}, { isHidden: true }),
};
export type BodyLabelColorThemeParams = typeof BodyLabelColorThemeParams

export function BodyLabelColorTheme(ctx: ThemeDataContext, props: PD.Values<BodyLabelColorThemeParams>): ColorTheme<BodyLabelColorThemeParams> {
    const volume = ctx.volume;
    const store = volume && BodyLabels.get(volume);

    if (!volume || !store) {
        return {
            factory: BodyLabelColorTheme,
            granularity: 'uniform',
            color: () => props.unassignedColor,
            props,
            description: Description,
        };
    }

    const colors = new Array<Color>(MaxBodyId + 1).fill(props.unassignedColor);
    for (const body of store.bodies) colors[body.id] = body.color;
    const labelAt = makeLabelAtPosition(volume, store.labels);

    const color: LocationColor = (location: Location) => {
        if (!isPositionLocation(location)) return props.unassignedColor;
        return colors[labelAt(location.position)];
    };

    return {
        factory: BodyLabelColorTheme,
        granularity: 'vertex',
        preferSmoothing: false,
        color,
        props,
        description: Description,
        contextHash: store.version,
        legend: TableLegend(store.bodies.map(b => [b.name, b.color] as [string, Color])),
    };
}

export const BodyLabelColorThemeProvider: ColorTheme.Provider<BodyLabelColorThemeParams, 'body-label'> = {
    name: 'body-label',
    label: 'Body Label',
    category: ColorThemeCategory.Misc,
    factory: BodyLabelColorTheme,
    getParams: () => BodyLabelColorThemeParams,
    defaultValues: PD.getDefaultValues(BodyLabelColorThemeParams),
    // Applicable to any volume so the theme survives param normalisation before labels exist.
    isApplicable: (ctx: ThemeDataContext) => !!ctx.volume,
};
