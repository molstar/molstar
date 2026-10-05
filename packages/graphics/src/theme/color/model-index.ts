/**
 * Copyright (c) 2022-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jason Pattle <jpattle@exscientia.co.uk>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '@molstar/core/util/color';
import type { Location } from '@molstar/model/model/location';
import { StructureElement, Bond, Model } from '@molstar/model/model/structure';
import type { ColorTheme, LocationColor } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { getPaletteParams, getPalette } from '@molstar/core/util/color/palette';
import type { TableLegend, ScaleLegend } from '@molstar/core/util/legend';
import { ColorThemeCategory } from './categories.js';

const DefaultColor = Color(0xcccccc);
const Description = 'Gives every model a unique color based on its index.';

export const ModelIndexColorThemeParams = {
  ...getPaletteParams({ type: 'colors', colorList: 'many-distinct' }),
};
export type ModelIndexColorThemeParams = typeof ModelIndexColorThemeParams;
export function getModelIndexColorThemeParams(ctx: ThemeDataContext) {
  return PD.clone(ModelIndexColorThemeParams);
}

export function ModelIndexColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<ModelIndexColorThemeParams>,
): ColorTheme<ModelIndexColorThemeParams> {
  let color: LocationColor;
  let legend: ScaleLegend | TableLegend | undefined;
  let contextHash = -1;

  if (ctx.structure) {
    // max-index is the same for all models
    const size = (Model.MaxIndex.get(ctx.structure.models[0]).value ?? -1) + 1;
    contextHash = size;

    const palette = getPalette(size, props);
    legend = palette.legend;

    color = (location: Location): Color => {
      if (StructureElement.Location.is(location)) {
        return palette.color(Model.Index.get(location.unit.model).value || 0)!;
      } else if (Bond.isLocation(location)) {
        return palette.color(Model.Index.get(location.aUnit.model).value || 0)!;
      }
      return DefaultColor;
    };
  } else {
    color = () => DefaultColor;
  }

  return {
    factory: ModelIndexColorTheme,
    granularity: 'instance',
    color,
    props,
    contextHash,
    description: Description,
    legend,
  };
}

export const ModelIndexColorThemeProvider: ColorTheme.Provider<ModelIndexColorThemeParams, 'model-index'> = {
  name: 'model-index',
  label: 'Model Index',
  category: ColorThemeCategory.Chain,
  factory: ModelIndexColorTheme,
  getParams: getModelIndexColorThemeParams,
  defaultValues: PD.getDefaultValues(ModelIndexColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.elementCount > 0,
};
