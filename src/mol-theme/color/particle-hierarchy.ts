/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import type { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';
import { DistinctColorsParams, getDistinctBaseColors, getDistinctGroupColors } from '../../mol-util/color/distinct';
import { getPalette } from '../../mol-util/color/palette';
import { ColorThemeCategory } from './categories';
import { Particle, ParticleList } from '../../mol-model/particles/particle-list';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every particle compartment a distinct base color and every entity within a compartment a distinct variation of that color.';

export const ParticleHierarchyColorThemeParams = {
    by: PD.Select('name', PD.arrayToOptions(['name', 'function'] as const), { description: 'Give a distinct color to every entity name or to every entity function' }),
    variability: PD.Numeric(20, { min: 1, max: 180, step: 1 }, { description: 'Hue range of the entity colors around the compartment color' }),
    shift: PD.Numeric(0, { min: 0, max: 100, step: 1 }, { description: 'Rotate the entity colors within a compartment' }),
    compartmentPalette: PD.MappedStatic('auto', {
        auto: PD.Group({
            compartmentShift: PD.Numeric(0, { min: 0, max: 100, step: 1 }, { description: 'Rotate the compartment base colors' }),
        }, { isFlat: true }),
        colors: PD.Group({
            list: PD.ColorList('set-1'),
        }, { isFlat: true }),
        generate: PD.Group({
            ...DistinctColorsParams,
            maxCount: PD.Numeric(75, { min: 1, max: 250, step: 1 }),
        }, { isFlat: true }),
    }, {
        options: [
            ['auto', 'Auto'],
            ['colors', 'Color List'],
            ['generate', 'Generate Distinct'],
        ],
        description: 'Base colors of the compartments, or of the entities if there are no compartments'
    }),
};
export type ParticleHierarchyColorThemeParams = typeof ParticleHierarchyColorThemeParams
type Props = PD.Values<ParticleHierarchyColorThemeParams>

function getParticleList(ctx: ThemeDataContext): ParticleList | undefined {
    if (ctx.particles) return ctx.particles;
    return undefined;
}

interface HierarchyColors {
    /** Per-particle index into `colors`, -1 for particles without a compartment. */
    colorIndex: Int32Array
    colors: Color[]
    /** Base color and label per compartment (per entity group if there are no compartments). */
    baseColors: Color[]
    baseLabels: string[]
}

interface EntityGroups {
    /** Per-particle entity group index, -1 for particles without a group. */
    groups: ArrayLike<number>
    getLabel: (group: number) => string
}

function getBaseColors(count: number, props: Props): Color[] {
    const p = props.compartmentPalette;
    if (p.name === 'auto') return getDistinctBaseColors(count, p.params.compartmentShift);

    const { color } = getPalette(count, { palette: p });
    const colors: Color[] = [];
    for (let i = 0; i < count; ++i) colors.push(color(i));
    return colors;
}

function getEntityGroups(particles: ParticleList, by: Props['by']): EntityGroups | undefined {
    const { count, entities, entityInfo } = particles;
    if (!entities) return undefined;

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
        return { groups, getLabel: g => functionNames[g] ?? `Function ${g}` };
    }

    return { groups: entities, getLabel: g => entityInfo?.get(g)?.name ?? `Entity ${g}` };
}

function buildEntityColors(particles: ParticleList, entityGroups: EntityGroups | undefined, props: Props): HierarchyColors {
    const { count } = particles;
    const colorIndex = new Int32Array(count).fill(-1);
    const labels: string[] = [];
    if (!entityGroups) return { colorIndex, colors: [], baseColors: [], baseLabels: labels };

    const { groups, getLabel } = entityGroups;
    const groupMap = new Map<number, number>();
    for (let i = 0; i < count; ++i) {
        const g = groups[i];
        if (g >= 0 && !groupMap.has(g)) {
            groupMap.set(g, groupMap.size);
            labels.push(getLabel(g));
        }
    }
    for (let i = 0; i < count; ++i) {
        const g = groups[i];
        if (g >= 0) colorIndex[i] = groupMap.get(g)!;
    }

    const colors = getBaseColors(groupMap.size, props);
    return { colorIndex, colors, baseColors: colors, baseLabels: labels };
}

function buildHierarchyColors(particles: ParticleList, props: Props): HierarchyColors {
    const { count, compartments, compartmentInfo } = particles;
    const entityGroups = getEntityGroups(particles, props.by);
    if (!compartments) return buildEntityColors(particles, entityGroups, props);

    const groups = entityGroups?.groups;
    const colorIndex = new Int32Array(count).fill(-1);
    const colors: Color[] = [];
    const baseLabels: string[] = [];

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

        const e = groups ? groups[i] : -1;
        if (e < 0) {
            hasNoEntity[ci] = true;
        } else if (!entityMaps[ci].has(e)) {
            entityMaps[ci].set(e, entityMaps[ci].size);
        }
    }

    const compartmentCount = compartmentMap.size;
    const baseColors = getBaseColors(compartmentCount, props);
    const entityOffsets = new Int32Array(compartmentCount);
    const noEntitySlots = new Int32Array(compartmentCount).fill(-1);

    for (let ci = 0; ci < compartmentCount; ++ci) {
        const compartmentLabel = compartmentInfo?.get(compartmentIds[ci])?.name ?? `Compartment ${compartmentIds[ci]}`;
        baseLabels.push(compartmentLabel);

        if (hasNoEntity[ci]) {
            noEntitySlots[ci] = colors.length;
            colors.push(baseColors[ci]);
        }

        entityOffsets[ci] = colors.length;
        const entityCount = entityMaps[ci].size;
        const groupColors = getDistinctGroupColors(entityCount, baseColors[ci], props.variability, props.shift);
        for (let ei = 0; ei < entityCount; ++ei) {
            colors.push(groupColors[ei]);
        }
    }

    for (let i = 0; i < count; ++i) {
        const c = compartments[i];
        if (c < 0) continue;

        const ci = compartmentMap.get(c)!;
        const e = groups ? groups[i] : -1;
        colorIndex[i] = e < 0
            ? noEntitySlots[ci]
            : entityOffsets[ci] + entityMaps[ci].get(e)!;
    }

    return { colorIndex, colors, baseColors, baseLabels };
}

export function getParticleHierarchyColorThemeParams(ctx: ThemeDataContext) {
    return PD.clone(ParticleHierarchyColorThemeParams);
}

export function ParticleHierarchyColorTheme(ctx: ThemeDataContext, props: Props): ColorTheme<ParticleHierarchyColorThemeParams> {
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
        legend
    };
}

export const ParticleHierarchyColorThemeProvider: ColorTheme.Provider<ParticleHierarchyColorThemeParams, 'particle-hierarchy'> = {
    name: 'particle-hierarchy',
    label: 'Particle Hierarchy',
    category: ColorThemeCategory.Particle,
    factory: ParticleHierarchyColorTheme,
    getParams: getParticleHierarchyColorThemeParams,
    defaultValues: PD.getDefaultValues(ParticleHierarchyColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!getParticleList(ctx)
};
