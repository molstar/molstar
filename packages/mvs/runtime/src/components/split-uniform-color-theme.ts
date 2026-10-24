/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ColorTheme } from '@molstar/graphics/theme/color';
import { ColorThemeCategory } from '@molstar/graphics/theme/color/categories';
import type { ThemeDataContext } from '@molstar/graphics/theme/theme';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { SplitColorProp } from '@molstar/core/util/color/value-palette-params';
import { isMVSStructure } from './is-mvs-model-prop.js';

export { type OptionalSplitColorProp, SplitColorProp } from '@molstar/core/util/color/value-palette-params';

export const MVSSplitUniformColorThemeParams = {
  value: SplitColorProp.PD(),
};
export type MVSSplitUniformColorThemeParams = typeof MVSSplitUniformColorThemeParams;
export type MVSSplitUniformColorThemeProps = PD.Values<MVSSplitUniformColorThemeParams>;

export function MVSSplitUniformColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<MVSSplitUniformColorThemeParams>,
): ColorTheme<MVSSplitUniformColorThemeParams> {
  const [color1, color2] = SplitColorProp.toTuple(props.value);

  return {
    factory: MVSSplitUniformColorTheme,
    granularity: 'groupInstance',
    color: (location, isSecondary) => (isSecondary ? color2 : color1),
    props: props,
    description:
      'Gives everything the same, uniform color, optionally split into two colors (when representation supports split coloring, e.g. for carbohydrate symbols, nucleic acid cartoon)',
  };
}

export const MVSSplitUniformColorThemeProvider: ColorTheme.Provider<
  MVSSplitUniformColorThemeParams,
  'mvs-split-uniform'
> = {
  name: 'mvs-split-uniform',
  label: 'MVS Split Uniform',
  category: ColorThemeCategory.Misc,
  factory: MVSSplitUniformColorTheme,
  getParams: () => MVSSplitUniformColorThemeParams,
  defaultValues: PD.getDefaultValues(MVSSplitUniformColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && isMVSStructure(ctx.structure),
};
