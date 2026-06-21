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
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';
import { ColorLists, getColorListFromName } from '../../mol-util/color/lists';
import { ColorThemeCategory } from './categories';
import { Particle, ParticleList } from '../../mol-model/particles/particle-list';
import { Structure, StructureElement, Bond } from '../../mol-model/structure';
import { Volume } from '../../mol-model/volume';

const DefaultList = 'dark-2';
const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every unique particle entity a distinct color.';

export const ParticleEntityColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type ParticleEntityColorThemeParams = typeof ParticleEntityColorThemeParams

function getParticleList(ctx: ThemeDataContext): ParticleList | undefined {
    if (ctx.particles) return ctx.particles;
    if (ctx.structure) {
        const p = Structure.ParticleList.get(ctx.structure);
        if (p) return p;
    }
    if (ctx.volume) {
        const p = Volume.ParticleList.get(ctx.volume);
        if (p) return p;
    }
    return undefined;
}

function buildEntityColorIndex(particles: ParticleList): { colorIndex: Int32Array, entityCount: number } {
    const { count, entities } = particles;
    const colorIndex = new Int32Array(count).fill(-1);
    if (!entities) return { colorIndex, entityCount: 0 };

    // Map raw entity indices to dense sequential palette indices.
    const entitySet = new Map<number, number>();
    for (let i = 0; i < count; ++i) {
        const e = entities[i];
        if (e >= 0 && !entitySet.has(e)) entitySet.set(e, entitySet.size);
    }
    for (let i = 0; i < count; ++i) {
        const e = entities[i];
        if (e >= 0) colorIndex[i] = entitySet.get(e)!;
    }
    return { colorIndex, entityCount: entitySet.size };
}

export function getParticleEntityColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(ParticleEntityColorThemeParams);
    const particles = getParticleList(ctx);
    if (particles) {
        const { entityCount } = buildEntityColorIndex(particles);
        if (entityCount > ColorLists[DefaultList].list.length) {
            params.palette.defaultValue.name = 'colors';
            params.palette.defaultValue.params = {
                ...params.palette.defaultValue.params,
                list: { kind: 'interpolate', colors: getColorListFromName(DefaultList).list }
            };
        }
    }
    return params;
}

export function ParticleEntityColorTheme(ctx: ThemeDataContext, props: PD.Values<ParticleEntityColorThemeParams>): ColorTheme<ParticleEntityColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    const particles = getParticleList(ctx);
    if (particles) {
        const { colorIndex, entityCount } = buildEntityColorIndex(particles);
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
            if (StructureElement.Location.is(location)) {
                return pick(location.unit.conformation.operator.group);
            }
            if (Bond.isLocation(location)) {
                return pick(location.aUnit.conformation.operator.group);
            }
            if (Volume.Cell.isLocation(location)) {
                const inst = location.volume.instances[location.instance];
                return pick(inst?.group ?? -1);
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
        legend
    };
}

export const ParticleEntityColorThemeProvider: ColorTheme.Provider<ParticleEntityColorThemeParams, 'particle-entity'> = {
    name: 'particle-entity',
    label: 'Particle Entity',
    category: ColorThemeCategory.Particle,
    factory: ParticleEntityColorTheme,
    getParams: getParticleEntityColorThemeParams,
    defaultValues: PD.getDefaultValues(ParticleEntityColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!getParticleList(ctx)?.entities
};
