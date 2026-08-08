/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorScale } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import type { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { ScaleLegend, TableLegend } from '../../mol-util/legend';
import { ColorThemeCategory } from './categories';
import { Particle, ParticleList } from '../../mol-model/particles/particle-list';
import { Structure, StructureElement, Bond } from '../../mol-model/structure';
import { Volume } from '../../mol-model/volume';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Colors particles by a per-particle scalar attribute (e.g. score, class, CC).';

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

export const ParticleAttributeColorThemeParams = {
    attribute: PD.Select('', [['', 'None']] as [string, string][]),
    domain: PD.MappedStatic('auto', {
        auto: PD.Group({}),
        custom: PD.Interval([0, 1], { step: 0.001 }),
    }),
    list: PD.ColorList('red-white-blue', { presetKind: 'scale' }),
};
export type ParticleAttributeColorThemeParams = typeof ParticleAttributeColorThemeParams

export function getParticleAttributeColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(ParticleAttributeColorThemeParams);
    const particles = getParticleList(ctx);
    if (particles?.attributeInfo && particles.attributeInfo.size > 0) {
        const options: [string, string][] = [];
        particles.attributeInfo.forEach((info, key) => options.push([key, info.label]));
        params.attribute = PD.Select(options[0][0], options);
    }
    return params;
}

export function ParticleAttributeColorTheme(ctx: ThemeDataContext, props: PD.Values<ParticleAttributeColorThemeParams>): ColorTheme<ParticleAttributeColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    const particles = getParticleList(ctx);
    if (particles?.attributes && particles.attributeInfo) {
        const attrData = particles.attributes.get(props.attribute);
        const attrInfo = particles.attributeInfo.get(props.attribute);

        if (attrData && attrInfo) {
            const domain: [number, number] = props.domain.name === 'custom'
                ? props.domain.params as [number, number]
                : [attrInfo.min, attrInfo.max];

            const scale = ColorScale.create({
                reverse: false,
                domain,
                listOrName: props.list.colors,
            });
            legend = scale.legend;

            const pick = (index: number): Color => {
                if (index < 0 || index >= attrData.length) return DefaultColor;
                return scale.color(attrData[index]);
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
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: ParticleAttributeColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend,
    };
}

export const ParticleAttributeColorThemeProvider: ColorTheme.Provider<ParticleAttributeColorThemeParams, 'particle-attribute'> = {
    name: 'particle-attribute',
    label: 'Particle Attribute',
    category: ColorThemeCategory.Particle,
    factory: ParticleAttributeColorTheme,
    getParams: getParticleAttributeColorThemeParams,
    defaultValues: PD.getDefaultValues(ParticleAttributeColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => {
        const p = getParticleList(ctx);
        return !!p?.attributes && p.attributes.size > 0;
    },
};
