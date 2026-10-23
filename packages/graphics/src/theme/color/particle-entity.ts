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
import { Particle, type ParticleList } from '@molstar/model/model/particles/particle-list';

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
  if (!entities) return 0;

  if (by === 'function') {
    const functions = new Set<string>();
    entityInfo?.forEach((info) => {
      if (info.function) functions.add(info.function);
    });
    return functions.size;
  }

  if (entityInfo) return entityInfo.size;

  const entitySet = new Set<number>();
  for (let i = 0, il = entities.length; i < il; ++i) {
    if (entities[i] >= 0) entitySet.add(entities[i]);
  }
  return entitySet.size;
}

interface EntityColorIndex {
  colorIndex: Int32Array;
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
  const colorIndex = new Int32Array(count).fill(-1);
  if (!entities) return { colorIndex, entityCount: 0 };

  if (by === 'function') {
    // Entities sharing a function get the same color, entities without one stay uncolored.
    const functionSet = new Map<string, number>();
    const entityToColorIdx = new Map<number, number>();
    entityInfo?.forEach((info, entityIdx) => {
      if (entityIdx < 0 || info.function === undefined) return;
      if (!functionSet.has(info.function)) functionSet.set(info.function, functionSet.size);
      entityToColorIdx.set(entityIdx, functionSet.get(info.function)!);
    });
    for (let i = 0; i < count; ++i) {
      const ci = entityToColorIdx.get(entities[i]);
      if (ci !== undefined) colorIndex[i] = ci;
    }
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
    for (let i = 0; i < count; ++i) {
      const e = entities[i];
      if (e >= 0 && !entitySet.has(e)) entitySet.set(e, entitySet.size);
    }
  }
  for (let i = 0; i < count; ++i) {
    const e = entities[i];
    if (e >= 0) colorIndex[i] = entitySet.get(e) ?? -1;
  }
  return { colorIndex, entityCount: entitySet.size };
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

  const particles = getParticleList(ctx);
  if (particles) {
    const { colorIndex, entityCount } = getEntityColorIndex(particles, props.by);
    const palette = getPalette(entityCount, props);
    legend = palette.legend;

    const pick = (particleIndex: number) => {
      const ci = colorIndex[particleIndex];
      return ci >= 0 ? palette.color(ci) : DefaultColor;
    };

    color = (location: Location): Color => {
      if (Particle.isLocation(location)) {
        return pick(location.index);
      }
      return DefaultColor;
    };
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
