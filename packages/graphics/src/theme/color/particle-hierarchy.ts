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
import { TableLegend, type ScaleLegend } from '@molstar/core/util/legend';
import { DistinctColorsParams, getDistinctBaseColors, getDistinctGroupColors } from '@molstar/core/util/color/distinct';
import { getPalette } from '@molstar/core/util/color/palette';
import { ColorThemeCategory } from './categories.js';
import { Particle, type ParticleList } from '@molstar/model/model/particles/particle-list';

const DefaultColor = Color(0xcccccc);
const Description =
  'Gives every particle compartment a distinct base color and every entity within a compartment a distinct variation of that color.';

export const ParticleHierarchyColorThemeParams = {
  by: PD.Select('name', PD.arrayToOptions(['name', 'function'] as const), {
    description: 'Give a distinct color to every entity name or to every entity function',
  }),
  variability: PD.Numeric(
    20,
    { min: 1, max: 180, step: 1 },
    { description: 'Hue range of the entity colors around the compartment color' },
  ),
  shift: PD.Numeric(0, { min: 0, max: 100, step: 1 }, { description: 'Rotate the entity colors within a compartment' }),
  compartmentPalette: PD.MappedStatic(
    'auto',
    {
      auto: PD.Group(
        {
          compartmentShift: PD.Numeric(
            0,
            { min: 0, max: 100, step: 1 },
            { description: 'Rotate the compartment base colors' },
          ),
        },
        { isFlat: true },
      ),
      colors: PD.Group(
        {
          list: PD.ColorList('set-1'),
        },
        { isFlat: true },
      ),
      generate: PD.Group(
        {
          ...DistinctColorsParams,
          maxCount: PD.Numeric(75, { min: 1, max: 250, step: 1 }),
        },
        { isFlat: true },
      ),
    },
    {
      options: [
        ['auto', 'Auto'],
        ['colors', 'Color List'],
        ['generate', 'Generate Distinct'],
      ],
      description: 'Base colors of the compartments',
    },
  ),
};
export type ParticleHierarchyColorThemeParams = typeof ParticleHierarchyColorThemeParams;
type Props = PD.Values<ParticleHierarchyColorThemeParams>;

type HierarchicalParticleList = ParticleList & Required<Pick<ParticleList, 'compartments' | 'entities'>>;

function getParticleList(ctx: ThemeDataContext): HierarchicalParticleList | undefined {
  if (ctx.particles?.compartments && ctx.particles.entities) return ctx.particles as HierarchicalParticleList;
  return undefined;
}

interface HierarchyColors {
  /** Per-particle index into `colors`, -1 for particles without a compartment. */
  colorIndex: Int32Array;
  colors: Color[];
  /** Base color and label per compartment. */
  baseColors: Color[];
  baseLabels: string[];
}

interface EntityGroups {
  /** Per-particle entity group index, -1 for particles without a group. */
  groups: ArrayLike<number>;
  getLabel: (group: number) => string;
}

function getBaseColors(count: number, props: Props): Color[] {
  const p = props.compartmentPalette;
  if (p.name === 'auto') return getDistinctBaseColors(count, p.params.compartmentShift);

  const { color } = getPalette(count, { palette: p });
  const colors: Color[] = [];
  for (let i = 0; i < count; ++i) colors.push(color(i));
  return colors;
}

function getEntityGroups(particles: HierarchicalParticleList, by: Props['by']): EntityGroups {
  const { count, entities, entityInfo } = particles;

  if (by === 'function') {
    // Entities sharing a function are grouped together, entities without one stay ungrouped.
    const functionNames: string[] = [];
    const functionIndex = new Map<string, number>();
    const entityToGroup = new Map<number, number>();
    entityInfo?.forEach((info, entityIdx) => {
      if (entityIdx < 0 || info.function === undefined) return;
      let gi = functionIndex.get(info.function);
      if (gi === undefined) {
        gi = functionNames.length;
        functionIndex.set(info.function, gi);
        functionNames.push(info.function);
      }
      entityToGroup.set(entityIdx, gi);
    });

    const groups = new Int32Array(count).fill(-1);
    for (let i = 0; i < count; ++i) {
      const gi = entityToGroup.get(entities[i]);
      if (gi !== undefined) groups[i] = gi;
    }
    return { groups, getLabel: (g) => functionNames[g] ?? `Function ${g}` };
  }

  return { groups: entities, getLabel: (g) => entityInfo?.get(g)?.name ?? `Entity ${g}` };
}

/** What the colors index by: per particle, a slot in the color list; per compartment, its slots. */
interface HierarchyIndex {
  /** Per-particle index into the color list, -1 for particles without a compartment. */
  colorIndex: Int32Array;
  /** Per dense compartment: its raw index, whether it has particles without an entity, its distinct entity count. */
  compartmentIds: number[];
  hasNoEntity: boolean[];
  entityCounts: number[];
}

/**
 * The O(particles) part, the same for every palette: per particle list and `by`, made once. The particle transforms
 * create the theme on every update (perf-audit-2 09/P6). The arrays of a list do not change; its `entityInfo` may
 * grow (Simularium frames share one and add types met later), so the index is remade when its size changed.
 */
const hierarchyIndexCache = new WeakMap<
  ParticleList,
  Map<Props['by'], { entityInfoSize: number; index: HierarchyIndex }>
>();

function getHierarchyIndex(particles: HierarchicalParticleList, by: Props['by']): HierarchyIndex {
  let byList = hierarchyIndexCache.get(particles);
  if (!byList) hierarchyIndexCache.set(particles, (byList = new Map()));
  const entityInfoSize = particles.entityInfo?.size ?? 0;
  let cached = byList.get(by);
  if (!cached || cached.entityInfoSize !== entityInfoSize) {
    byList.set(by, (cached = { entityInfoSize, index: buildHierarchyIndex(particles, by) }));
  }
  return cached.index;
}

function buildHierarchyIndex(particles: HierarchicalParticleList, by: Props['by']): HierarchyIndex {
  const { count, compartments } = particles;
  const groups = getEntityGroups(particles, by).groups;
  const colorIndex = new Int32Array(count).fill(-1);

  // Map raw compartment indices to dense sequential indices and collect the
  // distinct entities of each compartment.
  const compartmentMap = new Map<number, number>();
  const compartmentIds: number[] = [];
  const entityMaps: Map<number, number>[] = [];
  const hasNoEntity: boolean[] = [];
  for (let i = 0; i < count; ++i) {
    const c = compartments[i];
    if (c < 0) continue;

    let ci = compartmentMap.get(c);
    if (ci === undefined) {
      ci = compartmentMap.size;
      compartmentMap.set(c, ci);
      compartmentIds.push(c);
      entityMaps.push(new Map());
      hasNoEntity.push(false);
    }

    const e = groups[i];
    if (e < 0) {
      hasNoEntity[ci] = true;
    } else if (!entityMaps[ci].has(e)) {
      entityMaps[ci].set(e, entityMaps[ci].size);
    }
  }

  // the color list: per compartment its base color (if it has particles without an entity), then its entities
  const compartmentCount = compartmentMap.size;
  const entityOffsets = new Int32Array(compartmentCount);
  const noEntitySlots = new Int32Array(compartmentCount).fill(-1);
  const entityCounts: number[] = [];
  let slots = 0;
  for (let ci = 0; ci < compartmentCount; ++ci) {
    if (hasNoEntity[ci]) noEntitySlots[ci] = slots++;
    entityOffsets[ci] = slots;
    entityCounts.push(entityMaps[ci].size);
    slots += entityMaps[ci].size;
  }

  for (let i = 0; i < count; ++i) {
    const c = compartments[i];
    if (c < 0) continue;

    const ci = compartmentMap.get(c)!;
    const e = groups[i];
    colorIndex[i] = e < 0 ? noEntitySlots[ci] : entityOffsets[ci] + entityMaps[ci].get(e)!;
  }

  return { colorIndex, compartmentIds, hasNoEntity, entityCounts };
}

function buildHierarchyColors(particles: HierarchicalParticleList, props: Props): HierarchyColors {
  const { compartmentInfo } = particles;
  const { colorIndex, compartmentIds, hasNoEntity, entityCounts } = getHierarchyIndex(particles, props.by);

  const colors: Color[] = [];
  const baseLabels: string[] = [];
  const compartmentCount = compartmentIds.length;
  const baseColors = getBaseColors(compartmentCount, props);

  for (let ci = 0; ci < compartmentCount; ++ci) {
    const compartmentLabel = compartmentInfo?.get(compartmentIds[ci])?.name ?? `Compartment ${compartmentIds[ci]}`;
    baseLabels.push(compartmentLabel);

    if (hasNoEntity[ci]) colors.push(baseColors[ci]);

    const entityCount = entityCounts[ci];
    const groupColors = getDistinctGroupColors(entityCount, baseColors[ci], props.variability, props.shift);
    for (let ei = 0; ei < entityCount; ++ei) {
      colors.push(groupColors[ei]);
    }
  }

  return { colorIndex, colors, baseColors, baseLabels };
}

export function getParticleHierarchyColorThemeParams(ctx: ThemeDataContext) {
  return ParticleHierarchyColorThemeParams;
}

export function ParticleHierarchyColorTheme(
  ctx: ThemeDataContext,
  props: Props,
): ColorTheme<ParticleHierarchyColorThemeParams> {
  let color: LocationColor;
  let legend: ScaleLegend | TableLegend | undefined;

  const particles = getParticleList(ctx);
  if (particles) {
    const { colorIndex, colors, baseColors, baseLabels } = buildHierarchyColors(particles, props);

    legend = TableLegend(baseColors.map((c, i) => [baseLabels[i], c] as [string, Color]));

    color = (location: Location): Color => {
      if (Particle.isLocation(location)) {
        const ci = colorIndex[location.index];
        return ci >= 0 ? colors[ci] : DefaultColor;
      }
      return DefaultColor;
    };
  } else {
    color = () => DefaultColor;
  }

  return {
    factory: ParticleHierarchyColorTheme,
    granularity: 'instance',
    color,
    props,
    description: Description,
    legend,
  };
}

export const ParticleHierarchyColorThemeProvider: ColorTheme.Provider<
  ParticleHierarchyColorThemeParams,
  'particle-hierarchy'
> = {
  name: 'particle-hierarchy',
  label: 'Particle Hierarchy',
  category: ColorThemeCategory.Particle,
  factory: ParticleHierarchyColorTheme,
  getParams: getParticleHierarchyColorThemeParams,
  defaultValues: PD.getDefaultValues(ParticleHierarchyColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!getParticleList(ctx),
};
