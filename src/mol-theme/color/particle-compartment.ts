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
const Description = 'Gives every unique particle compartment a distinct color.';

export const ParticleCompartmentColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type ParticleCompartmentColorThemeParams = typeof ParticleCompartmentColorThemeParams

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

function buildCompartmentColorIndex(particles: ParticleList): { colorIndex: Int32Array, compartmentCount: number } {
    const { count, compartments } = particles;
    const colorIndex = new Int32Array(count).fill(-1);
    if (!compartments) return { colorIndex, compartmentCount: 0 };

    // Map raw compartment indices to dense sequential palette indices.
    const compartmentSet = new Map<number, number>();
    for (let i = 0; i < count; ++i) {
        const c = compartments[i];
        if (c >= 0 && !compartmentSet.has(c)) compartmentSet.set(c, compartmentSet.size);
    }
    for (let i = 0; i < count; ++i) {
        const c = compartments[i];
        if (c >= 0) colorIndex[i] = compartmentSet.get(c)!;
    }
    return { colorIndex, compartmentCount: compartmentSet.size };
}

export function getParticleCompartmentColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(ParticleCompartmentColorThemeParams);
    const particles = getParticleList(ctx);
    if (particles) {
        const { compartmentCount } = buildCompartmentColorIndex(particles);
        if (compartmentCount > ColorLists[DefaultList].list.length) {
            params.palette.defaultValue.name = 'colors';
            params.palette.defaultValue.params = {
                ...params.palette.defaultValue.params,
                list: { kind: 'interpolate', colors: getColorListFromName(DefaultList).list }
            };
        }
    }
    return params;
}

export function ParticleCompartmentColorTheme(ctx: ThemeDataContext, props: PD.Values<ParticleCompartmentColorThemeParams>): ColorTheme<ParticleCompartmentColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    const particles = getParticleList(ctx);
    if (particles) {
        const { colorIndex, compartmentCount } = buildCompartmentColorIndex(particles);
        const palette = getPalette(compartmentCount, props);
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
        factory: ParticleCompartmentColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    };
}

export const ParticleCompartmentColorThemeProvider: ColorTheme.Provider<ParticleCompartmentColorThemeParams, 'particle-compartment'> = {
    name: 'particle-compartment',
    label: 'Particle Compartment',
    category: ColorThemeCategory.Particle,
    factory: ParticleCompartmentColorTheme,
    getParams: getParticleCompartmentColorThemeParams,
    defaultValues: PD.getDefaultValues(ParticleCompartmentColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!getParticleList(ctx)
};
