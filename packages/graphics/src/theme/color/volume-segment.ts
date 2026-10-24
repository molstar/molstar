/**
 * Copyright (c) 2022-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
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

const DefaultColor = Color(0xcccccc);
const Description = 'Gives every volume segment a unique color.';

export const VolumeSegmentColorThemeParams = {
  ...getPaletteParams({ type: 'colors', colorList: 'many-distinct' }),
};
export type VolumeSegmentColorThemeParams = typeof VolumeSegmentColorThemeParams;
export function getVolumeSegmentColorThemeParams(ctx: ThemeDataContext) {
  return PD.clone(VolumeSegmentColorThemeParams);
}

export function VolumeSegmentColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<VolumeSegmentColorThemeParams>,
): ColorTheme<VolumeSegmentColorThemeParams> {
  let color: LocationColor;
  let legend: ScaleLegend | TableLegend | undefined;

  const segmentation = ctx.volume && Volume.Segmentation.get(ctx.volume);

  if (segmentation) {
    const size = segmentation.segments.size;
    // a segment's place among the segments, its palette index (a Map: indexOf was O(segments) per location, perf-audit-2
    // 09/P18)
    const places = new Map<number, number>();
    for (const segment of segmentation.segments.keys()) places.set(segment, places.size);

    const palette = getPalette(size, props);
    legend = palette.legend;

    color = (location: Location): Color => {
      if (Volume.Segment.isLocation(location)) {
        return palette.color(places.get(location.segment) ?? -1);
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
    legend,
  };
}

export const VolumeSegmentColorThemeProvider: ColorTheme.Provider<VolumeSegmentColorThemeParams, 'volume-segment'> = {
  name: 'volume-segment',
  label: 'Volume Segment',
  category: ColorThemeCategory.Misc,
  factory: VolumeSegmentColorTheme,
  getParams: getVolumeSegmentColorThemeParams,
  defaultValues: PD.getDefaultValues(VolumeSegmentColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.volume && !!Volume.Segmentation.get(ctx.volume),
};
