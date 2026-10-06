/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Location } from '@molstar/model/model/location';
import { Bond, StructureElement, Unit } from '@molstar/model/model/structure';
import type { ColorTheme, LocationColor } from '@molstar/graphics/theme/color';
import { ColorThemeCategory } from '@molstar/graphics/theme/color/categories';
import type { ThemeDataContext } from '@molstar/graphics/theme/theme';
import { Color } from '@molstar/core/util/color';
import { getPalette, getPaletteParams } from '@molstar/core/util/color/palette';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import { SIFTSMapping } from '@molstar/model/props/sequence/sifts-mapping';

const DefaultColor = Color(0xfafafa);
const Description = 'Assigns a color based on SIFTS mapping.';

// same colors for same accessions
const globalAccessionMap = new Map<string, number>();

export const SIFTSMappingColorThemeParams = {
  ...getPaletteParams({ type: 'colors', colorList: 'set-1' }),
};
export type SIFTSMappingColorThemeParams = typeof SIFTSMappingColorThemeParams;
export function getSIFTSMappingColorThemeParams(ctx: ThemeDataContext) {
  return SIFTSMappingColorThemeParams; // TODO return copy
}
export function SIFTSMappingColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<SIFTSMappingColorThemeParams>,
): ColorTheme<SIFTSMappingColorThemeParams> {
  let color: LocationColor;

  if (ctx.structure) {
    for (const m of ctx.structure.models) {
      const mapping = SIFTSMapping.Provider.get(m).value;
      if (!mapping) continue;
      for (const acc of mapping.accession) {
        if (!acc || globalAccessionMap.has(acc)) continue;
        globalAccessionMap.set(acc, globalAccessionMap.size);
      }
    }

    const l = StructureElement.Location.create(ctx.structure);
    const palette = getPalette(globalAccessionMap.size + 1, props, { valueLabel: (i) => `${i}` });
    const colorMap = new Map<string, Color>();

    const getColor = (location: StructureElement.Location) => {
      const key = SIFTSMapping.getKey(location);
      if (!key) return DefaultColor;

      if (colorMap.has(key)) return colorMap.get(key)!;

      const color = palette.color(globalAccessionMap.get(key)!);
      colorMap.set(key, color);
      return color;
    };

    color = (location: Location): Color => {
      if (StructureElement.Location.is(location) && Unit.isAtomic(location.unit)) {
        return getColor(location);
      } else if (Bond.isLocation(location)) {
        l.unit = location.aUnit;
        l.element = location.aUnit.elements[location.aIndex];
        return getColor(l);
      }
      return DefaultColor;
    };
  } else {
    color = () => DefaultColor;
  }

  return {
    factory: SIFTSMappingColorTheme,
    granularity: 'group',
    preferSmoothing: true,
    color,
    props,
    description: Description,
  };
}

export const SIFTSMappingColorThemeProvider: ColorTheme.Provider<SIFTSMappingColorThemeParams, 'sifts-mapping'> = {
  name: 'sifts-mapping',
  label: 'SIFTS Mapping',
  category: ColorThemeCategory.Residue,
  factory: SIFTSMappingColorTheme,
  getParams: getSIFTSMappingColorThemeParams,
  defaultValues: PD.getDefaultValues(SIFTSMappingColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure?.models.some((m) => SIFTSMapping.Provider.isApplicable(m)),
  ensureCustomProperties: {
    attach: async (ctx: CustomProperty.Context, data: ThemeDataContext) => {
      if (!data.structure) return;

      for (const m of data.structure.models) {
        await SIFTSMapping.Provider.attach(ctx, m, void 0, true);
      }
    },
    detach: (data) => {
      if (!data.structure) return;
      for (const m of data.structure.models) {
        SIFTSMapping.Provider.ref(m, false);
      }
    },
  },
};
