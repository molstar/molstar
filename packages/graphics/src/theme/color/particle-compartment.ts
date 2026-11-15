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
import { fillByIndex, fillByLocation, forEachPageList, getPagedPainting, slotPagesHave } from './particle-pages.js';

import {
  Particle,
  type ParticleList,
  type ParticleRuns,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';

const DefaultList = 'dark-2';
const DefaultColor = Color(0xcccccc);
const Description = 'Gives every unique particle compartment a distinct color.';

export const ParticleCompartmentColorThemeParams = {
  ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type ParticleCompartmentColorThemeParams = typeof ParticleCompartmentColorThemeParams;

function getParticleList(ctx: ThemeDataContext): ParticleList | undefined {
  if (ctx.particles) return ctx.particles;
  return undefined;
}

/**
 * A particle's palette index: its compartment's, dense in the order of `compartmentInfo` (a slot list's and a view's of
 * some tiles the same: the catalog's), those it does not hold after them as first met, -1 without one; in the
 * compartments' form (`ParticleValues`: runs a run, not a particle)
 */
function buildCompartmentColorIndex(particles: ParticleList): {
  colorIndex: ParticleValues<Int32Array> | undefined;
  compartmentCount: number;
} {
  const { count, compartments, compartmentInfo } = particles;
  if (!compartments) return { colorIndex: undefined, compartmentCount: 0 };
  const compartmentSet = new Map<number, number>();
  if (compartmentInfo) {
    for (const c of compartmentInfo.keys())
      if (c >= 0 && !compartmentSet.has(c)) compartmentSet.set(c, compartmentSet.size);
  }
  ParticleValues.forEachRun(compartments, count, (_start, _end, c) => {
    if (c >= 0 && !compartmentSet.has(c)) compartmentSet.set(c, compartmentSet.size);
  });
  const colorIndex = ParticleValues.map(compartments, count, (c) => (c >= 0 ? compartmentSet.get(c)! : -1));
  return { colorIndex, compartmentCount: compartmentSet.size };
}

export function getParticleCompartmentColorThemeParams(ctx: ThemeDataContext) {
  const params = PD.clone(ParticleCompartmentColorThemeParams);
  const particles = getParticleList(ctx);
  if (particles) {
    const compartmentCount = particles.compartmentInfo?.size ?? 0;
    if (compartmentCount > ColorLists[DefaultList].list.length) {
      params.palette.defaultValue.name = 'colors';
      params.palette.defaultValue.params = {
        ...params.palette.defaultValue.params,
        list: { kind: 'interpolate', colors: getColorListFromName(DefaultList).list },
      };
    }
  }
  return params;
}

export function ParticleCompartmentColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<ParticleCompartmentColorThemeParams>,
): ColorTheme<ParticleCompartmentColorThemeParams> {
  let color: LocationColor;
  let legend: ScaleLegend | TableLegend | undefined;
  let particleRuns: ParticleRuns | undefined;
  let fillParticleColors: ColorTheme<ParticleCompartmentColorThemeParams>['fillParticleColors'];

  const particles = getParticleList(ctx);
  if (particles?.pages) {
    // a slot list: the compartments by its info (one for all pages), else by those its pages hold, in their order
    const dense = new Map<number, number>();
    if (particles.compartmentInfo) {
      for (const c of particles.compartmentInfo.keys()) if (c >= 0 && !dense.has(c)) dense.set(c, dense.size);
    } else {
      forEachPageList(particles, (list) => {
        if (!list.compartments) return;
        ParticleValues.forEachRun(list.compartments, list.count, (_s, _e, c) => {
          if (c >= 0 && !dense.has(c)) dense.set(c, dense.size);
        });
      });
    }
    const palette = getPalette(dense.size, props);
    legend = palette.legend;
    const colors: Color[] = [];
    for (let k = 0; k < dense.size; ++k) colors[k] = palette.color(k);
    const painting = getPagedPainting((page) => {
      const { compartments } = page;
      const colorIndex = compartments
        ? ParticleValues.map(compartments, page.count, (c) => (c >= 0 ? (dense.get(c) ?? -1) : -1))
        : undefined;
      return {
        color: (i) => {
          const ci = colorIndex ? ParticleValues.get(colorIndex, i) : -1;
          return ci >= 0 ? colors[ci] : DefaultColor;
        },
        fill: (out, offset) => fillByIndex(colorIndex, colors, page.count, DefaultColor, out, offset),
      };
    }, DefaultColor);
    color = painting.color;
    fillParticleColors = painting.fill;
  } else if (particles) {
    const { colorIndex, compartmentCount } = buildCompartmentColorIndex(particles);
    // a color a run of the index (a pack's: a run a target)
    if (colorIndex && ParticleValues.isRuns(colorIndex)) particleRuns = colorIndex;
    const palette = getPalette(compartmentCount, props);
    legend = palette.legend;

    // a color per compartment, not per particle (the palette interpolates)
    const colors: Color[] = [];
    for (let k = 0; k < compartmentCount; ++k) colors[k] = palette.color(k);
    const pick = (particleIndex: number) => {
      const ci = colorIndex ? ParticleValues.get(colorIndex, particleIndex) : -1;
      return ci >= 0 ? colors[ci] : DefaultColor;
    };

    color = (location: Location): Color => {
      if (Particle.isLocation(location)) {
        return pick(location.index);
      }
      return DefaultColor;
    };
    fillParticleColors = (list, out, offset) =>
      list === particles
        ? fillByIndex(colorIndex, colors, list.count, DefaultColor, out, offset)
        : fillByLocation(color, list, out, offset);
  } else {
    color = () => DefaultColor;
  }

  return {
    factory: ParticleCompartmentColorTheme,
    granularity: 'instance',
    color,
    props,
    description: Description,
    legend,
    particleRuns,
    fillParticleColors,
  };
}

export const ParticleCompartmentColorThemeProvider: ColorTheme.Provider<
  ParticleCompartmentColorThemeParams,
  'particle-compartment'
> = {
  name: 'particle-compartment',
  label: 'Particle Compartment',
  category: ColorThemeCategory.Particle,
  factory: ParticleCompartmentColorTheme,
  getParams: getParticleCompartmentColorThemeParams,
  defaultValues: PD.getDefaultValues(ParticleCompartmentColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) =>
    !!getParticleList(ctx) &&
    (ctx.particles?.compartments !== undefined || slotPagesHave(ctx.particles!, 'compartments')),
};
