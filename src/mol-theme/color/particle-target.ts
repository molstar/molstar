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
const Description = 'Gives every unique particle target a distinct color.';

export const ParticleTargetColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type ParticleTargetColorThemeParams = typeof ParticleTargetColorThemeParams

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

function buildTargetColorIndex(particles: ParticleList): { colorIndex: Int32Array, targetCount: number } {
    const { count, targets } = particles;
    const targetSet = new Map<number, number>();
    for (let i = 0; i < count; ++i) {
        const t = targets[i];
        if (!targetSet.has(t)) targetSet.set(t, targetSet.size);
    }
    const colorIndex = new Int32Array(count);
    for (let i = 0; i < count; ++i) {
        colorIndex[i] = targetSet.get(targets[i])!;
    }
    return { colorIndex, targetCount: targetSet.size };
}

export function getParticleTargetColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(ParticleTargetColorThemeParams);
    const particles = getParticleList(ctx);
    if (particles) {
        const { targetCount } = buildTargetColorIndex(particles);
        if (targetCount > ColorLists[DefaultList].list.length) {
            params.palette.defaultValue.name = 'colors';
            params.palette.defaultValue.params = {
                ...params.palette.defaultValue.params,
                list: { kind: 'interpolate', colors: getColorListFromName(DefaultList).list }
            };
        }
    }
    return params;
}

export function ParticleTargetColorTheme(ctx: ThemeDataContext, props: PD.Values<ParticleTargetColorThemeParams>): ColorTheme<ParticleTargetColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    const particles = getParticleList(ctx);
    if (particles) {
        const { colorIndex, targetCount } = buildTargetColorIndex(particles);
        const palette = getPalette(targetCount, props);
        legend = palette.legend;

        const pick = (particleIndex: number) => {
            const ci = colorIndex[particleIndex];
            return ci !== undefined ? palette.color(ci) : DefaultColor;
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
        factory: ParticleTargetColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    };
}

export const ParticleTargetColorThemeProvider: ColorTheme.Provider<ParticleTargetColorThemeParams, 'particle-target'> = {
    name: 'particle-target',
    label: 'Particle Target',
    category: ColorThemeCategory.Particle,
    factory: ParticleTargetColorTheme,
    getParams: getParticleTargetColorThemeParams,
    defaultValues: PD.getDefaultValues(ParticleTargetColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!getParticleList(ctx)
};
