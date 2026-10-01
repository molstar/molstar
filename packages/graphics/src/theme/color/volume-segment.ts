/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '@molstar/core/util/color';
import type { Location } from '@molstar/model/model/location';
import type { ColorTheme, LocationColor } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { getPaletteParams, getPalette } from '@molstar/core/util/color/palette';
import type { TableLegend, ScaleLegend } from '@molstar/core/util/legend';
import { Volume } from '@molstar/model/model/volume/volume';
import { ColorThemeCategory } from './categories.js';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every volume segment a unique color.';

export const VolumeSegmentColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: 'many-distinct' }),
};
export type VolumeSegmentColorThemeParams = typeof VolumeSegmentColorThemeParams
export function getVolumeSegmentColorThemeParams(ctx: ThemeDataContext) {
    return PD.clone(VolumeSegmentColorThemeParams);
}

export function VolumeSegmentColorTheme(ctx: ThemeDataContext, props: PD.Values<VolumeSegmentColorThemeParams>): ColorTheme<VolumeSegmentColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    const segmentation = ctx.volume && Volume.Segmentation.get(ctx.volume);

    if (segmentation) {
        const size = segmentation.segments.size;
        const segments = Array.from(segmentation.segments.keys());

        const palette = getPalette(size, props);
        legend = palette.legend;

        color = (location: Location): Color => {
            if (Volume.Segment.isLocation(location)) {
                return palette.color(segments.indexOf(location.segment));
            }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: VolumeSegmentColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    };
}

export const VolumeSegmentColorThemeProvider: ColorTheme.Provider<VolumeSegmentColorThemeParams, 'volume-segment'> = {
    name: 'volume-segment',
    label: 'Volume Segment',
    category: ColorThemeCategory.Misc,
    factory: VolumeSegmentColorTheme,
    getParams: getVolumeSegmentColorThemeParams,
    defaultValues: PD.getDefaultValues(VolumeSegmentColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.volume && !!Volume.Segmentation.get(ctx.volume)
};