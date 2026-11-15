/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { ColorLists, getColorListFromName } from '@molstar/core/util/color/lists';
import { ColorThemeCategory } from './categories.js';
import { fillByLocation, setParticleRgba } from './particle-pages.js';
import { Particle, type ParticleList } from '@molstar/model/model/particles/particle-list';
import { ParticleTileView } from '@molstar/model/model/particles/particle-tile-view';

const DefaultList = 'dark-2';
const DefaultColor = Color(0xcccccc);
const Description = 'Gives every particle a unique color based on its index in the particle list.';

/**
 * A particle's index and the count it is of: in a view of tiles its index in the whole (its key) of the set's count, so
 * that it keeps its color in views of other tiles; else its index in the list
 */
/** the interpolated palette's steps a table holds: far finer than a channel's 256 levels over its few colors */
const LutSize = 65536;
const lutCache = new Map<string, Color[]>();

/**
 * An interpolated palette's colors of `count` indices as a table of `LutSize` steps where there are more (a slot list's
 * millions: a palette's interpolation a particle about 20 ns), made once for the palette and count; none for fewer, or
 * a palette of cycled colors
 */
function getInterpolatedTable(
  count: number,
  props: PD.Values<ParticleIndexColorThemeParams>,
  color: (i: number) => Color,
) {
  const p = props.palette;
  if (count <= LutSize || p.name !== 'colors' || p.params.list.kind !== 'interpolate') return undefined;
  const key = `${count}|${JSON.stringify(p)}`;
  let lut = lutCache.get(key);
  if (!lut) {
    lut = new Array(LutSize);
    for (let k = 0; k < LutSize; ++k) lut[k] = color(Math.round((k * (count - 1)) / (LutSize - 1)));
    if (lutCache.size > 16) lutCache.clear();
    lutCache.set(key, lut);
  }
  return lut;
}

/** `lut`'s colors as opaque RGBA words (the default color last), made once a table */
const rgbaTables = new WeakMap<Color[], Uint32Array>();
function getRgbaTable(lut: Color[]) {
  let table = rgbaTables.get(lut);
  if (!table) {
    table = new Uint32Array(lut.length + 1);
    const bytes = new Uint8Array(table.buffer);
    for (let k = 0; k <= lut.length; ++k) setParticleRgba(bytes, k * 4, k < lut.length ? lut[k] : DefaultColor);
    rgbaTables.set(lut, table);
  }
  return table;
}

function getIndexing(particles: ParticleList) {
  // a slot list: a particle's slot (a tile view's its key), of the slots
  if (particles.pages) return { count: particles.count, indexOf: (i: number) => i };
  const view = ParticleTileView.of(particles);
  const keys = particles.keys;
  if (view && keys) return { count: view.set.count, indexOf: (i: number) => keys[i] };
  return { count: particles.count, indexOf: (i: number) => i };
}

export const ParticleIndexColorThemeParams = {
  ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type ParticleIndexColorThemeParams = typeof ParticleIndexColorThemeParams;

function getParticleList(ctx: ThemeDataContext): ParticleList | undefined {
  if (ctx.particles) return ctx.particles;
  return undefined;
}

export function getParticleIndexColorThemeParams(ctx: ThemeDataContext) {
  const params = PD.clone(ParticleIndexColorThemeParams);
  const particles = getParticleList(ctx);
  if (particles) {
    const { count } = getIndexing(particles);
    if (count > ColorLists[DefaultList].list.length) {
      params.palette.defaultValue.name = 'colors';
      params.palette.defaultValue.params = {
        ...params.palette.defaultValue.params,
        list: { kind: 'interpolate', colors: getColorListFromName(DefaultList).list },
      };
    }
  }
  return params;
}

export function ParticleIndexColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<ParticleIndexColorThemeParams>,
): ColorTheme<ParticleIndexColorThemeParams> {
  let color: LocationColor;
  let legend: ScaleLegend | TableLegend | undefined;
  let fillParticleColors: ColorTheme<ParticleIndexColorThemeParams>['fillParticleColors'];

  const particles = getParticleList(ctx);
  if (particles) {
    const { count, indexOf } = getIndexing(particles);
    const palette = getPalette(count, props);
    legend = palette.legend;

    const lut = getInterpolatedTable(count, props, palette.color);
    // a table's step by an index scaled once (the fill's the same)
    const scale = (LutSize - 1) / (count - 1);
    const pick = lut
      ? (index: number) => (index >= 0 && index < count ? lut[(index * scale + 0.5) | 0] : DefaultColor)
      : (index: number) => (index >= 0 && index < count ? palette.color(index) : DefaultColor);

    color = (location: Location): Color => {
      if (Particle.isLocation(location)) {
        return pick(indexOf(location.index));
      }
      return DefaultColor;
    };
    // a list's colors by its particles' indices in the whole: a slot list's page's (its keys), the list's own
    fillParticleColors = (list, out, offset) => {
      const keys = list.keys;
      if (list === particles && !particles.pages) {
        for (let i = 0; i < list.count; ++i) setParticleRgba(out, offset + i * 4, pick(indexOf(i)));
      } else if (keys && (particles.pages || ParticleTileView.of(particles))) {
        if (lut && offset % 4 === 0 && out.byteOffset % 4 === 0) {
          // the table's colors as RGBA words, a particle's by its key scaled once (a tile's block: a word a particle)
          const words = new Uint32Array(out.buffer, out.byteOffset + offset, list.count);
          const table = getRgbaTable(lut);
          for (let i = 0; i < list.count; ++i) {
            const k = keys[i];
            words[i] = k >= 0 && k < count ? table[(k * scale + 0.5) | 0] : table[LutSize];
          }
        } else {
          for (let i = 0; i < list.count; ++i) setParticleRgba(out, offset + i * 4, pick(keys[i]));
        }
      } else {
        fillByLocation(color, list, out, offset);
      }
    };
  } else {
    color = () => DefaultColor;
  }

  return {
    factory: ParticleIndexColorTheme,
    granularity: 'instance',
    color,
    props,
    description: Description,
    legend,
    fillParticleColors,
  };
}

export const ParticleIndexColorThemeProvider: ColorTheme.Provider<ParticleIndexColorThemeParams, 'particle-index'> = {
  name: 'particle-index',
  label: 'Particle Index',
  category: ColorThemeCategory.Particle,
  factory: ParticleIndexColorTheme,
  getParams: getParticleIndexColorThemeParams,
  defaultValues: PD.getDefaultValues(ParticleIndexColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!getParticleList(ctx),
};
