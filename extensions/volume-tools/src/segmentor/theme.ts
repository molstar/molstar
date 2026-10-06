/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { ColorTheme, type LocationColor } from '@molstar/graphics/theme/color';
import { ColorThemeCategory } from '@molstar/graphics/theme/color/categories';
import type { ThemeDataContext } from '@molstar/graphics/theme/theme';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Color } from '@molstar/core/util/color';
import { TableLegend } from '@molstar/core/util/legend';
import type { Location } from '@molstar/model/model/location';
import { isPositionLocation } from '@molstar/graphics/geo/util/location-iterator';
import { makeLabelAtPosition } from '@molstar/volume-tools-extension/voxel-labels';
import { BodyLabels } from './labels.js';
import { MaxBodyId } from './types.js';

const Description = 'Colors a volume surface by the body each voxel is assigned to.';

export const BodyLabelColorThemeParams = {
  unassignedColor: PD.Color(Color(0x9a9a9a), { description: 'Color of voxels not assigned to any body.' }),
  /** Mirrors `LabelStore.version`; bumping it forces a color update after labels change. */
  version: PD.Numeric(0, {}, { isHidden: true }),
};
export type BodyLabelColorThemeParams = typeof BodyLabelColorThemeParams;

export function BodyLabelColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<BodyLabelColorThemeParams>,
): ColorTheme<BodyLabelColorThemeParams> {
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
    legend: TableLegend(store.bodies.map((b) => [b.name, b.color] as [string, Color])),
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
