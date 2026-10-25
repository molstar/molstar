/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import type { Location } from '@molstar/model/model/location';
import type {
  ColorType,
  ColorTypeDirect,
  ColorTypeGrid,
  ColorTypeLocation,
} from '@molstar/graphics/geo/geometry/color-data';
import { deepEqual } from '@molstar/core/util';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type ThemeDataContext, ThemeRegistry, type ThemeProvider } from './theme.js';
import type { BuiltInColorThemes } from './color/catalog.js';
import type { ScaleLegend, TableLegend } from '@molstar/core/util/legend';
import type { Texture, TextureFilter } from '@molstar/graphics/gl/texture';
import type { Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { ColorThemeCategory } from './color/categories.js';
import type { ColorListEntry } from '@molstar/core/util/color/color';
import { getPrecision } from '@molstar/core/util/number';
import { SortedArray } from '@molstar/core/data/int/sorted-array';
import { normalize } from '@molstar/core/math/interpolate';

export type LocationColor = (location: Location, isSecondary: boolean) => Color;

export interface ColorVolume {
  colors: Texture;
  dimension: Vec3;
  transform: Vec4;
}

export { ColorTheme };

type ColorThemeShared<P extends PD.Params, G extends ColorType> = {
  readonly factory: ColorTheme.Factory<P, G>;
  readonly props: Readonly<PD.Values<P>>;
  /**
   * if palette is defined, 24bit RGB color value normalized to interval [0, 1]
   * is used as index to the colors
   */
  readonly palette?: Readonly<ColorTheme.Palette>;
  readonly preferSmoothing?: boolean;
  readonly contextHash?: number;
  readonly description?: string;
  readonly legend?: Readonly<ScaleLegend | TableLegend>;
};

type ColorThemeLocation<P extends PD.Params> = {
  readonly granularity: ColorTypeLocation;
  readonly color: LocationColor;
} & ColorThemeShared<P, ColorTypeLocation>;

type ColorThemeGrid<P extends PD.Params> = {
  readonly granularity: ColorTypeGrid;
  readonly grid: ColorVolume;
} & ColorThemeShared<P, ColorTypeGrid>;

type ColorThemeDirect<P extends PD.Params> = {
  readonly granularity: ColorTypeDirect;
} & ColorThemeShared<P, ColorTypeDirect>;

type ColorTheme<P extends PD.Params, G extends ColorType = ColorTypeLocation> = G extends ColorTypeLocation
  ? ColorThemeLocation<P>
  : G extends ColorTypeGrid
    ? ColorThemeGrid<P>
    : G extends ColorTypeDirect
      ? ColorThemeDirect<P>
      : never;

namespace ColorTheme {
  export const Category = ColorThemeCategory;

  export interface Palette {
    colors: Color[];
    filter?: TextureFilter;
    domain?: [number, number];
    defaultColor?: Color;
  }

  /** Samples of an interpolated palette with offsets; the most samples of one (the texture's width). */
  const PaletteInterpolatedSamples = 512;
  const MaxPaletteSamples = 4096;

  export function Palette(
    list: ColorListEntry[],
    kind: 'set' | 'interpolate',
    domain?: [number, number],
    defaultColor?: Color,
  ): Palette {
    const colors: Color[] = [];

    const hasOffsets = list.every((c) => Array.isArray(c));
    if (hasOffsets) {
      const sorted = [...list] as [Color, number][];
      sorted.sort((a, b) => a[1] - b[1]);
      const src = sorted.map((c) => c[0]);
      const values = SortedArray.ofSortedArray(sorted.map((c) => c[1]));

      if (kind === 'interpolate') {
        // Interpolated between the stops at their offsets (as `ColorScale` does), sampled finely enough for the
        // texture's linear filter; before the first stop and after the last, their colors.
        const count = PaletteInterpolatedSamples;
        for (let i = 0; i < count; ++i) {
          const t = normalize(i, 0, count - 1);
          const j = SortedArray.findPredecessorIndex(values, t);
          if (j <= 0) colors[i] = src[0];
          else if (j >= src.length) colors[i] = src[src.length - 1];
          else {
            const a = values[j - 1];
            const b = values[j];
            colors[i] = Color.interpolate(src[j - 1], src[j], b > a ? (t - a) / (b - a) : 1);
          }
        }
      } else {
        // A set: each stop's color up to the middle between it and the next. Offsets' precision sets the samples
        // (capped: an offset with many decimals does not make a huge texture; at least 2: integer offsets have no
        // decimals). After the last stop, its color (perf-audit-2 09/B3: `undefined`, black).
        let maxPrecision = 0;
        for (const e of list) {
          if (Array.isArray(e)) {
            const p = getPrecision(e[1]);
            if (p > maxPrecision) maxPrecision = p;
          }
        }
        const count = Math.max(2, Math.min(Math.pow(10, maxPrecision), MaxPaletteSamples));

        const _off: number[] = [];
        for (let i = 0, il = values.length - 1; i < il; ++i) {
          _off.push(values[i] + (values[i + 1] - values[i]) / 2);
        }
        _off.push(values[values.length - 1]);
        const off = SortedArray.ofSortedArray(_off);

        for (let i = 0, il = Math.max(count, list.length); i < il; ++i) {
          const t = normalize(i, 0, count - 1);
          const j = Math.min(SortedArray.findPredecessorIndex(off, t), src.length - 1);
          colors[i] = src[j];
        }
      }
    } else {
      for (const e of list) {
        if (Array.isArray(e)) colors.push(e[0]);
        else colors.push(e);
      }
    }

    return {
      colors,
      filter: kind === 'set' ? 'nearest' : 'linear',
      domain,
      defaultColor,
    };
  }

  export const PaletteScale = (1 << 24) - 2; // reserve (1 << 24) - 1 for undefiend values

  export type Props = { [k: string]: any };
  export type Factory<P extends PD.Params, G extends ColorType> = (
    ctx: ThemeDataContext,
    props: PD.Values<P>,
  ) => ColorTheme<P, G>;
  export const EmptyFactory = () => Empty;
  const EmptyColor = Color(0xcccccc);
  export const Empty: ColorTheme<{}> = {
    factory: EmptyFactory,
    granularity: 'uniform',
    color: () => EmptyColor,
    props: {},
  };

  /** per pair of themes whether they are equal: the visuals of a representation compare the same pair each update */
  const equalThemes = new WeakMap<ColorTheme<any, any>, WeakMap<ColorTheme<any, any>, boolean>>();

  export function areEqual(themeA: ColorTheme<any, any>, themeB: ColorTheme<any, any>) {
    if (themeA === themeB) return true;
    let byB = equalThemes.get(themeA);
    let equal = byB?.get(themeB);
    if (equal === undefined) {
      equal =
        themeA.contextHash === themeB.contextHash &&
        themeA.factory === themeB.factory &&
        deepEqual(themeA.props, themeB.props);
      if (!byB) equalThemes.set(themeA, (byB = new WeakMap()));
      byB.set(themeB, equal);
    }
    return equal;
  }

  export interface Provider<P extends PD.Params = any, Id extends string = string, G extends ColorType = ColorType>
    extends ThemeProvider<ColorTheme<P, G>, P, Id, G> {}
  export const EmptyProvider: Provider<{}> = {
    name: '',
    label: '',
    category: '',
    factory: EmptyFactory,
    getParams: () => ({}),
    defaultValues: {},
    isApplicable: () => true,
  };

  export type Registry = ThemeRegistry<ColorTheme<any, any>>;
  export function createRegistry() {
    const registry: Registry = new ThemeRegistry(EmptyProvider, 'color');
    return registry;
  }

  type _BuiltIn = typeof BuiltInColorThemes;
  export type BuiltIn = keyof _BuiltIn;
  export type ParamValues<C extends ColorTheme.Provider<any>> =
    C extends ColorTheme.Provider<infer P> ? PD.Values<P> : never;
  export type BuiltInParams<T extends BuiltIn> = Partial<ParamValues<_BuiltIn[T]>>;
}

export function ColorThemeProvider<P extends PD.Params, Id extends string>(
  p: ColorTheme.Provider<P, Id>,
): ColorTheme.Provider<P, Id> {
  return p;
}
