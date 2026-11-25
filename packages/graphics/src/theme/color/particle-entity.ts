/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { ColorLists, getColorListFromName } from '@molstar/core/util/color/lists';
import { ColorThemeCategory } from './categories.js';
import { fillByIndex, fillByLocation, getPagedPainting } from './particle-pages.js';
import {
  Particle,
  type ParticleList,
  type ParticleRuns,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';

const DefaultList = 'dark-2';
const DefaultColor = Color(0xcccccc);
const Description = 'Gives every unique particle entity a distinct color.';

export const ParticleEntityColorThemeParams = {
  by: PD.Select('name', PD.arrayToOptions(['name', 'function'] as const), {
    description: 'Give a distinct color to every entity name or to every entity function',
  }),
  ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type ParticleEntityColorThemeParams = typeof ParticleEntityColorThemeParams;
type Props = PD.Values<ParticleEntityColorThemeParams>;

function getParticleList(ctx: ThemeDataContext): ParticleList | undefined {
  if (ctx.particles) return ctx.particles;
  return undefined;
}

function getEntityCount(particles: ParticleList, by: Props['by']): number {
  const { entities, entityInfo } = particles;
  // a slot list's by its info (shared by its pages), its own entities none
  if (!entities && !particles.pages) return 0;

  if (by === 'function') {
    const functions = new Set<string>();
    entityInfo?.forEach((info) => {
      if (info.function) functions.add(info.function);
    });
    return functions.size;
  }

  if (entityInfo) return entityInfo.size;
  if (!entities) return 0;

  const entitySet = new Set<number>();
  ParticleValues.forEachRun(entities, particles.count, (_start, _end, e) => {
    if (e >= 0) entitySet.add(e);
  });
  return entitySet.size;
}

interface EntityColorIndex {
  /** a particle's palette index, -1 without; in the entities' form (`ParticleValues`: runs a run, not a particle) */
  colorIndex: ParticleValues<Int32Array> | undefined;
  entityCount: number;
}

/**
 * The O(particles) part, the same for every palette: per particle list and `by`, made once. The particle transforms
 * create the theme on every update (perf-audit-2 09/P6). The arrays of a list do not change; its `entityInfo` may
 * grow (Simularium frames share one and add types met later), so the index is remade when its size changed.
 */
const entityColorIndexCache = new WeakMap<
  ParticleList,
  Map<Props['by'], { entityInfoSize: number; index: EntityColorIndex }>
>();

function getEntityColorIndex(particles: ParticleList, by: Props['by']): EntityColorIndex {
  let byList = entityColorIndexCache.get(particles);
  if (!byList) entityColorIndexCache.set(particles, (byList = new Map()));
  const entityInfoSize = particles.entityInfo?.size ?? 0;
  let cached = byList.get(by);
  if (!cached || cached.entityInfoSize !== entityInfoSize) {
    byList.set(by, (cached = { entityInfoSize, index: buildEntityColorIndex(particles, by) }));
  }
  return cached.index;
}

function buildEntityColorIndex(particles: ParticleList, by: Props['by']): EntityColorIndex {
  const { count, entities, entityInfo } = particles;
  if (!entities) return { colorIndex: undefined, entityCount: 0 };

  if (by === 'function') {
    // Entities sharing a function get the same color, entities without one stay uncolored.
    const functionSet = new Map<string, number>();
    const entityToColorIdx = new Map<number, number>();
    entityInfo?.forEach((info, entityIdx) => {
      if (entityIdx < 0 || info.function === undefined) return;
      if (!functionSet.has(info.function)) functionSet.set(info.function, functionSet.size);
      entityToColorIdx.set(entityIdx, functionSet.get(info.function)!);
    });
    const colorIndex = mapByTable(entities, count, entityToColorIdx);
    return { colorIndex, entityCount: functionSet.size };
  }

  // Map raw entity indices to dense sequential palette indices.
  const entitySet = new Map<number, number>();
  if (entityInfo) {
    entityInfo.forEach((_, entityIdx) => {
      if (entityIdx >= 0 && !entitySet.has(entityIdx)) {
        entitySet.set(entityIdx, entitySet.size);
      }
    });
  } else {
    ParticleValues.forEachRun(entities, count, (_start, _end, e) => {
      if (e >= 0 && !entitySet.has(e)) entitySet.set(e, entitySet.size);
    });
  }
  const colorIndex = mapByTable(entities, count, entitySet);
  return { colorIndex, entityCount: entitySet.size };
}

/**
 * the entities' palette indices (-1 an entity not in `toIndex`), by a table of the entities' indices where they are
 * small (a map's lookup a run: the yeast tile-ordered 2.29 M runs)
 */
function mapByTable(entities: ParticleValues<Int32Array>, count: number, toIndex: Map<number, number>) {
  let max = -1;
  for (const e of toIndex.keys()) if (e > max) max = e;
  if (max >= 1 << 20) return ParticleValues.map(entities, count, (e) => (e >= 0 ? (toIndex.get(e) ?? -1) : -1));
  const table = new Int32Array(max + 1).fill(-1);
  for (const [e, ci] of toIndex) if (e >= 0) table[e] = ci;
  return ParticleValues.map(entities, count, (e) => (e >= 0 && e <= max ? table[e] : -1));
}

/**
 * A target's color by its entity (`ParticleTargetInfo.entity`: every particle of it has that entity), the palette's
 * index of the entity as the particles' (`buildEntityColorIndex`: by the entity info's order, or by function); none
 * without an entity info or a target's entity
 */
function getTargetColor(particles: ParticleList, by: Props['by'], colors: readonly Color[]) {
  const { entityInfo, targetInfo } = particles;
  if (!entityInfo) return undefined;
  const entityToColor = new Map<number, number>();
  if (by === 'function') {
    const functions = new Map<string, number>();
    entityInfo.forEach((info, e) => {
      if (e < 0 || info.function === undefined) return;
      if (!functions.has(info.function)) functions.set(info.function, functions.size);
      entityToColor.set(e, functions.get(info.function)!);
    });
  } else {
    entityInfo.forEach((_, e) => {
      if (e >= 0 && !entityToColor.has(e)) entityToColor.set(e, entityToColor.size);
    });
  }
  return (target: number): Color | undefined => {
    const e = targetInfo.get(target)?.entity;
    if (e === undefined) return undefined;
    const ci = entityToColor.get(e);
    return ci !== undefined && ci >= 0 && ci < colors.length ? colors[ci] : DefaultColor;
  };
}

export function getParticleEntityColorThemeParams(ctx: ThemeDataContext) {
  const params = PD.clone(ParticleEntityColorThemeParams);
  const particles = getParticleList(ctx);
  if (particles) {
    const entityCount = getEntityCount(particles, params.by.defaultValue);
    if (entityCount > ColorLists[DefaultList].list.length) {
      params.palette.defaultValue.name = 'colors';
      params.palette.defaultValue.params = {
        ...params.palette.defaultValue.params,
        list: { kind: 'interpolate', colors: getColorListFromName(DefaultList).list },
      };
    }
  }
  return params;
}

export function ParticleEntityColorTheme(
  ctx: ThemeDataContext,
  props: Props,
): ColorTheme<ParticleEntityColorThemeParams> {
  let color: LocationColor;
  let legend: ScaleLegend | TableLegend | undefined;
  let particleRuns: ParticleRuns | undefined;
  let fillParticleColors: ColorTheme<ParticleEntityColorThemeParams>['fillParticleColors'];
  let targetColor: ColorTheme<ParticleEntityColorThemeParams>['targetColor'];

  const particles = getParticleList(ctx);
  if (particles?.pages) {
    // a slot list: the palette by its info (its pages', one for all), a page's index its own (the same entities' order)
    const entityCount = getEntityCount(particles, props.by);
    const palette = getPalette(entityCount, props);
    legend = palette.legend;
    const colors: Color[] = [];
    for (let k = 0; k < entityCount; ++k) colors[k] = palette.color(k);
    const painting = getPagedPainting((page) => {
      const { colorIndex } = getEntityColorIndex(page, props.by);
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
    targetColor = getTargetColor(particles, props.by, colors);
  } else if (particles) {
    const { colorIndex, entityCount } = getEntityColorIndex(particles, props.by);
    // a color a run of the index (a pack's: a run a target)
    if (colorIndex && ParticleValues.isRuns(colorIndex)) particleRuns = colorIndex;
    const palette = getPalette(entityCount, props);
    legend = palette.legend;

    // a color per entity, not per particle (the palette interpolates)
    const colors: Color[] = [];
    for (let k = 0; k < entityCount; ++k) colors[k] = palette.color(k);
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
    if (particles.entityInfo) targetColor = getTargetColor(particles, props.by, colors);
  } else {
    color = () => DefaultColor;
  }

  return {
    factory: ParticleEntityColorTheme,
    granularity: 'instance',
    color,
    props,
    description: Description,
    legend,
    particleRuns,
    fillParticleColors,
    targetColor,
  };
}

export const ParticleEntityColorThemeProvider: ColorTheme.Provider<ParticleEntityColorThemeParams, 'particle-entity'> =
  {
    name: 'particle-entity',
    label: 'Particle Entity',
    category: ColorThemeCategory.Particle,
    factory: ParticleEntityColorTheme,
    getParams: getParticleEntityColorThemeParams,
    defaultValues: PD.getDefaultValues(ParticleEntityColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!getParticleList(ctx) && ctx.particles?.entityInfo !== undefined,
  };
