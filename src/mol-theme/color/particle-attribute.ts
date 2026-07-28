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

const DefaultColor = Color(0xCCCCCC);
const Description = 'Colors particles by a per-particle scalar attribute (e.g. score, class, CC).';

function getParticleList(ctx: ThemeDataContext): ParticleList | undefined {
    if (ctx.particles) return ctx.particles;
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
    if (particles?.attributes && particles.attributes.size > 0) {
        const options: [string, string][] = [];
        particles.attributes.forEach((attr, key) => {
            if (['number', 'float'].includes(attr.column.schema.valueType)) {
                options.push([key, attr.label]);
            }
        });
        if (options.length > 0) {
            params.attribute = PD.Select(options[0][0], options);
        }
    }
    return params;
}

export function ParticleAttributeColorTheme(ctx: ThemeDataContext, props: PD.Values<ParticleAttributeColorThemeParams>): ColorTheme<ParticleAttributeColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    const particles = getParticleList(ctx);
    if (particles?.attributes) {
        const attr = particles.attributes.get(props.attribute);

        if (attr) {
            const domain: [number, number] = props.domain.name === 'custom'
                ? props.domain.params as [number, number]
                : [attr.min, attr.max];

            const scale = ColorScale.create({
                reverse: false,
                domain,
                listOrName: props.list.colors,
            });
            legend = scale.legend;

            const pick = (index: number): Color => {
                if (index < 0 || index >= attr.column.rowCount) return DefaultColor;
                return scale.color(attr.column.value(index));
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
