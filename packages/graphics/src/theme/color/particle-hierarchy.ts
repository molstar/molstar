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
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { fillByIndex, fillByLocation, forEachPageList, getPagedPainting, slotPagesHave } from './particle-pages.js';

import {
  Particle,
  type ParticleList,
  type ParticleRuns,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';

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
  const particles = ctx.particles;
  if (particles?.compartments && particles.entities) return particles as HierarchicalParticleList;
  // a slot list: its pages' (each a flat list with both)
  if (particles && slotPagesHave(particles, 'compartments') && slotPagesHave(particles, 'entities')) {
    return particles as HierarchicalParticleList;
  }
  return undefined;
}

interface HierarchyColors {
  /** Per-particle index into `colors`, -1 for particles without a compartment (`ParticleValues`). */
  colorIndex: ParticleValues<Int32Array>;
  colors: Color[];
  /** Base color and label per compartment. */
  baseColors: Color[];
  baseLabels: string[];
}

interface EntityGroups {
  /** Per-particle entity group index, -1 for particles without a group (`ParticleValues`). */
  groups: ParticleValues<Int32Array>;
  getLabel: (group: number) => string;
}

/**
 * A compartment's entities' colors: deterministic and costly (distinct colors in Lab, about 100 ms for the yeast's
 * compartments), the theme made again every update: made once for its arguments
 */
const groupColorsCache = new Map<string, Color[]>();
function getGroupColors(count: number, base: Color, variability: number, shift: number): Color[] {
  const key = `${count}|${base}|${variability}|${shift}`;
  let colors = groupColorsCache.get(key);
  if (!colors) groupColorsCache.set(key, (colors = getDistinctGroupColors(count, base, variability, shift)));
  return colors;
}

function getBaseColors(count: number, props: Props): Color[] {
  const p = props.compartmentPalette;
  if (p.name === 'auto') return getDistinctBaseColors(count, p.params.compartmentShift);

  const { color } = getPalette(count, { palette: p });
  const colors: Color[] = [];
  for (let i = 0; i < count; ++i) colors.push(color(i));
  return colors;
}

/**
 * The groups of entities by function (`entityInfo`'s, in its order): an entity's group, the functions' names; entities
 * without one ungrouped
 */
function getFunctionGroups(entityInfo: ParticleList['entityInfo']) {
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
  return { functionNames, entityToGroup };
}

function getEntityGroups(particles: HierarchicalParticleList, by: Props['by']): EntityGroups {
  const { count, entities, entityInfo } = particles;

  if (by === 'function') {
    // Entities sharing a function are grouped together, entities without one stay ungrouped.
    const { functionNames, entityToGroup } = getFunctionGroups(entityInfo);
    const groups = ParticleValues.map(entities, count, (e) => entityToGroup.get(e) ?? -1);
    return { groups, getLabel: (g) => functionNames[g] ?? `Function ${g}` };
  }

  return { groups: entities, getLabel: (g) => entityInfo?.get(g)?.name ?? `Entity ${g}` };
}

/** What the colors index by: per particle, a slot in the color list; per compartment, its slots. */
interface HierarchyIndex {
  /** Per-particle index into the color list, -1 for particles without a compartment (`ParticleValues`). */
  colorIndex: ParticleValues<Int32Array>;
  /** Per dense compartment: its raw index, whether it has particles without an entity, its distinct entity count. */
  compartmentIds: number[];
  hasNoEntity: boolean[];
  entityCounts: number[];
}

/**
 * The O(particles) part, the same for every palette: per compartment array and `by`, made once (a list's copies share
 * it: one moved as a whole, `Particle.transform`, or with targets). The particle transforms create the theme on every
 * update (perf-audit-2 09/P6). The arrays of a list do not change; its `entityInfo` may grow (Simularium frames share
 * one and add types met later), so the index is remade when its size changed, or where the entities are others.
 */
const hierarchyIndexCache = new WeakMap<
  object,
  Map<
    Props['by'],
    {
      entities: ParticleValues<Int32Array>;
      entityInfo: unknown;
      entityInfoSize: number;
      compartmentInfo: unknown;
      index: HierarchyIndex;
    }
  >
>();

function getHierarchyIndex(particles: HierarchicalParticleList, by: Props['by']): HierarchyIndex {
  let byList = hierarchyIndexCache.get(particles.compartments);
  if (!byList) hierarchyIndexCache.set(particles.compartments, (byList = new Map()));
  const { entities, entityInfo, compartmentInfo } = particles;
  const entityInfoSize = entityInfo?.size ?? 0;
  let cached = byList.get(by);
  if (
    !cached ||
    cached.entityInfoSize !== entityInfoSize ||
    cached.entities !== entities ||
    cached.entityInfo !== entityInfo ||
    cached.compartmentInfo !== compartmentInfo
  ) {
    const index = buildHierarchyIndex(particles, by);
    byList.set(by, (cached = { entities, entityInfo, entityInfoSize, compartmentInfo, index }));
  }
  return cached.index;
}

/** The compartments and their entities of some lists (a slot list's pages: one for all), their colors' slots */
interface HierarchyStructure {
  compartmentMap: Map<number, number>;
  compartmentIds: number[];
  entityMaps: Map<number, number>[];
  hasNoEntity: boolean[];
  entityCounts: number[];
  noEntitySlots: Int32Array;
  entityOffsets: Int32Array;
}

/**
 * The compartments and their entities (or groups, `by`) in the catalog's order: the compartments by
 * `compartmentInfo`'s, an entity by `entityInfo`'s (a function by its first entity's), those the info does not hold
 * after them as met. Of the catalog where its compartments' info holds their entities (`ParticleCompartmentInfo.entities`,
 * a tile set's: a view of some tiles colors as the whole), else of `lists`' particles
 */
function collectHierarchy(
  lists: readonly HierarchicalParticleList[],
  by: Props['by'],
  catalog: ParticleList,
): HierarchyStructure {
  const { compartmentInfo, entityInfo } = catalog;
  const fromInfo =
    !!compartmentInfo && compartmentInfo.size > 0 && [...compartmentInfo.values()].every((c) => !!c.entities);
  if (!fromInfo) return makeHierarchy(lists, by, catalog, false);
  // of the catalog alone: made once a catalog (a following view's step creates the theme again: thousands of
  // entities sorted a step otherwise), again where its info is another or grew
  let byInfo = hierarchyStructureCache.get(compartmentInfo);
  if (!byInfo) hierarchyStructureCache.set(compartmentInfo, (byInfo = new Map()));
  const entityInfoSize = entityInfo?.size ?? 0;
  const compartmentInfoSize = compartmentInfo.size;
  const cached = byInfo.get(by);
  if (
    cached &&
    cached.entityInfo === entityInfo &&
    cached.entityInfoSize === entityInfoSize &&
    cached.compartmentInfoSize === compartmentInfoSize
  ) {
    return cached.structure;
  }
  const structure = makeHierarchy(lists, by, catalog, true);
  byInfo.set(by, { entityInfo, entityInfoSize, compartmentInfoSize, structure });
  return structure;
}

/** `collectHierarchy`'s of a catalog's info (its compartments holding their entities), by its `compartmentInfo` */
const hierarchyStructureCache = new WeakMap<
  object,
  Map<
    Props['by'],
    { entityInfo: unknown; entityInfoSize: number; compartmentInfoSize: number; structure: HierarchyStructure }
  >
>();

function makeHierarchy(
  lists: readonly HierarchicalParticleList[],
  by: Props['by'],
  catalog: ParticleList,
  fromInfo: boolean,
): HierarchyStructure {
  const { compartmentInfo, entityInfo } = catalog;
  /** per compartment (raw) its groups and whether it has particles without one, as met */
  const met = new Map<number, { groups: Set<number>; noEntity: boolean }>();
  const metOf = (c: number) => {
    let m = met.get(c);
    if (!m) met.set(c, (m = { groups: new Set(), noEntity: false }));
    return m;
  };
  if (fromInfo) {
    const toGroup = by === 'function' ? getFunctionGroups(entityInfo).entityToGroup : undefined;
    for (const [c, info] of compartmentInfo!) {
      if (c < 0) continue;
      const m = metOf(c);
      if (info.withoutEntity) m.noEntity = true;
      for (const e of info.entities!) {
        const g = toGroup ? (toGroup.get(e) ?? -1) : e;
        if (g < 0) m.noEntity = true;
        else m.groups.add(g);
      }
    }
  } else {
    for (const particles of lists) {
      const groups = getEntityGroups(particles, by).groups;
      // by the runs of both (`ParticleValues`): a run's lookups, not a particle's
      ParticleValues.forEachRun2(particles.compartments, groups, particles.count, (_start, _end, c, e) => {
        if (c < 0) return;
        const m = metOf(c);
        if (e < 0) m.noEntity = true;
        else m.groups.add(e);
      });
    }
  }

  // the order: a key's place in the info's, those it does not hold after them (as met: a stable sort)
  const rankIn = (keys: Iterable<number> | undefined) => {
    const rank = new Map<number, number>();
    if (keys) for (const k of keys) if (!rank.has(k)) rank.set(k, rank.size);
    return (k: number) => rank.get(k) ?? Number.MAX_SAFE_INTEGER;
  };
  const compartmentRank = rankIn(compartmentInfo?.keys());
  // a function group's id is its place among the functions (`getFunctionGroups`)
  const groupRank = by === 'function' ? (g: number) => g : rankIn(entityInfo?.keys());
  const compartmentIds = [...met.keys()].sort((a, b) => compartmentRank(a) - compartmentRank(b));
  const compartmentMap = new Map<number, number>();
  const entityMaps: Map<number, number>[] = [];
  const hasNoEntity: boolean[] = [];
  for (const c of compartmentIds) {
    const m = met.get(c)!;
    compartmentMap.set(c, compartmentMap.size);
    const sorted = [...m.groups].sort((a, b) => groupRank(a) - groupRank(b));
    entityMaps.push(new Map(sorted.map((g, k) => [g, k])));
    hasNoEntity.push(m.noEntity);
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
  return { compartmentMap, compartmentIds, entityMaps, hasNoEntity, entityCounts, noEntitySlots, entityOffsets };
}

/** `particles`' color index by `structure` (collected over it, or over the pages it is one of) */
function indexHierarchy(
  structure: HierarchyStructure,
  particles: HierarchicalParticleList,
  by: Props['by'],
): ParticleValues<Int32Array> {
  const { compartmentMap, entityMaps, noEntitySlots, entityOffsets } = structure;
  const { count, compartments } = particles;
  const groups = getEntityGroups(particles, by).groups;
  // a value a run, as runs where they are fewer than half the particles (a pack's lists: a run a target)
  const runEnds: number[] = [];
  const runValues: number[] = [];
  ParticleValues.forEachRun2(compartments, groups, count, (_start, end, c, e) => {
    const ci = c < 0 ? -1 : (compartmentMap.get(c) ?? -1);
    const slot = ci < 0 ? -1 : e < 0 ? noEntitySlots[ci] : entityOffsets[ci] + (entityMaps[ci].get(e) ?? -1);
    const value = ci >= 0 && e >= 0 && !entityMaps[ci].has(e) ? -1 : slot;
    if (runValues.length > 0 && runValues[runValues.length - 1] === value) runEnds[runEnds.length - 1] = end;
    else {
      runEnds.push(end);
      runValues.push(value);
    }
  });
  const offsets = new Int32Array(runEnds.length + 1);
  for (let r = 0; r < runEnds.length; ++r) offsets[r + 1] = runEnds[r];
  const runs = ParticleValues.ofRuns(offsets, Int32Array.from(runValues));
  return runValues.length * 2 <= count ? runs : ParticleValues.toArray(runs, count, Int32Array);
}

function buildHierarchyIndex(particles: HierarchicalParticleList, by: Props['by']): HierarchyIndex {
  const structure = collectHierarchy([particles], by, particles);
  const { compartmentIds, hasNoEntity, entityCounts } = structure;
  return { colorIndex: indexHierarchy(structure, particles, by), compartmentIds, hasNoEntity, entityCounts };
}

function buildHierarchyColors(particles: HierarchicalParticleList, props: Props): HierarchyColors {
  const { colorIndex, compartmentIds, hasNoEntity, entityCounts } = getHierarchyIndex(particles, props.by);
  return { colorIndex, ...buildColorList(particles, compartmentIds, hasNoEntity, entityCounts, props) };
}

/** The color list: per compartment its base color (where it has particles without an entity), then its entities' */
function buildColorList(
  particles: ParticleList,
  compartmentIds: number[],
  hasNoEntity: boolean[],
  entityCounts: number[],
  props: Props,
) {
  const { compartmentInfo } = particles;

  const colors: Color[] = [];
  const baseLabels: string[] = [];
  const compartmentCount = compartmentIds.length;
  const baseColors = getBaseColors(compartmentCount, props);

  for (let ci = 0; ci < compartmentCount; ++ci) {
    const compartmentLabel = compartmentInfo?.get(compartmentIds[ci])?.name ?? `Compartment ${compartmentIds[ci]}`;
    baseLabels.push(compartmentLabel);

    if (hasNoEntity[ci]) colors.push(baseColors[ci]);

    const entityCount = entityCounts[ci];
    const groupColors = getGroupColors(entityCount, baseColors[ci], props.variability, props.shift);
    for (let ei = 0; ei < entityCount; ++ei) {
      colors.push(groupColors[ei]);
    }
  }

  return { colors, baseColors, baseLabels };
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
  let particleRuns: ParticleRuns | undefined;
  let contextHash: number | undefined;
  let fillParticleColors: ColorTheme<ParticleHierarchyColorThemeParams>['fillParticleColors'];
  let targetColor: ColorTheme<ParticleHierarchyColorThemeParams>['targetColor'];

  const particles = getParticleList(ctx);
  if (particles?.pages) {
    // a slot list: the compartments and their entities in the catalog's order (`collectHierarchy`), a page's index by
    // them
    const lists: HierarchicalParticleList[] = [];
    forEachPageList(particles, (list) => {
      if (list.compartments && list.entities) lists.push(list as HierarchicalParticleList);
    });
    const structure = collectHierarchy(lists, props.by, particles);
    const { compartmentIds, hasNoEntity, entityCounts } = structure;
    // its colors those of the structure: another one (a tile bringing a compartment or an entity) colors all again
    contextHash = hashFnv32a([...compartmentIds, -1, ...entityCounts, -1, ...hasNoEntity.map((b) => (b ? 1 : 0))]);
    const { colors, baseColors, baseLabels } = buildColorList(
      particles,
      compartmentIds,
      hasNoEntity,
      entityCounts,
      props,
    );
    legend = TableLegend(baseColors.map((c, i) => [baseLabels[i], c] as [string, Color]));
    const painting = getPagedPainting((page) => {
      const colorIndex =
        page.compartments && page.entities
          ? indexHierarchy(structure, page as HierarchicalParticleList, props.by)
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
    targetColor = getTargetColor(particles, structure, colors, props.by);
  } else if (particles) {
    const { colorIndex, colors, baseColors, baseLabels } = buildHierarchyColors(particles, props);
    // a color a run of the index (a pack's: a target's particles one compartment and entity)
    if (ParticleValues.isRuns(colorIndex)) particleRuns = colorIndex;

    legend = TableLegend(baseColors.map((c, i) => [baseLabels[i], c] as [string, Color]));

    color = (location: Location): Color => {
      if (Particle.isLocation(location)) {
        const ci = ParticleValues.get(colorIndex, location.index);
        return ci >= 0 ? colors[ci] : DefaultColor;
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
    factory: ParticleHierarchyColorTheme,
    granularity: 'instance',
    color,
    props,
    description: Description,
    legend,
    particleRuns,
    fillParticleColors,
    targetColor,
    contextHash,
  };
}

/**
 * A target's color where its info gives its compartment (`ParticleTargetInfo.compartment`: its particles' all, a
 * pack's) and so its entity's: as `indexHierarchy` colors its particles; none otherwise (colored a particle)
 */
function getTargetColor(
  particles: ParticleList,
  structure: HierarchyStructure,
  colors: readonly Color[],
  by: Props['by'],
) {
  const { compartmentMap, entityMaps, noEntitySlots, entityOffsets } = structure;
  const entityToGroup = by === 'function' ? getFunctionGroups(particles.entityInfo).entityToGroup : undefined;
  return (target: number): Color | undefined => {
    const info = particles.targetInfo.get(target);
    const c = info?.compartment;
    if (c === undefined) return undefined;
    const ci = c < 0 ? -1 : (compartmentMap.get(c) ?? -1);
    if (ci < 0) return DefaultColor;
    const entity = info!.entity;
    const e = entity === undefined ? -1 : entityToGroup ? (entityToGroup.get(entity) ?? -1) : entity;
    if (e < 0) {
      const slot = noEntitySlots[ci];
      return slot >= 0 ? colors[slot] : DefaultColor;
    }
    const k = entityMaps[ci].get(e);
    return k === undefined ? DefaultColor : (colors[entityOffsets[ci] + k] ?? DefaultColor);
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
