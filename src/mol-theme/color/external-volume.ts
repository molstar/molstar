/**
 * Copyright (c) 2024-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Cai Huiyu <szmun.caihy@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorScale } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { Grid, Volume } from '../../mol-model/volume';
import { type PluginContext } from '../../mol-plugin/context';
import { isPositionLocation } from '../../mol-geo/util/location-iterator';
import { Vec3 } from '../../mol-math/linear-algebra';
import { clamp } from '../../mol-math/interpolate';
import { ColorThemeCategory } from './categories';

const Description = `Assigns a color based on volume value at a given vertex.`;

export const ExternalVolumeColorThemeParams = {
    volume: PD.ValueRef<Volume>(
        (ctx: PluginContext) => {
            const volumes = ctx.state.data.selectQ(q => q.root.subtree().filter(c => Volume.is(c.obj?.data)));
            return volumes.map(v => [v.transform.ref, v.obj?.label ?? '<unknown>'] as [string, string]);
        },
        (ref, getData) => getData(ref),
    ),
    coloring: PD.MappedStatic('absolute-value', {
        'absolute-value': PD.Group({
            domain: PD.MappedStatic('auto', {
                custom: PD.Interval([-1, 1], { step: 0.001 }),
                auto: PD.Group({
                    symmetric: PD.Boolean(false, { description: 'If true the automatic range is determined as [-|max|, |max|].' })
                })
            }),
            list: PD.ColorList('red-white-blue', { presetKind: 'scale' })
        }),
        'relative-value': PD.Group({
            domain: PD.MappedStatic('auto', {
                custom: PD.Interval([-1, 1], { step: 0.001 }),
                auto: PD.Group({
                    symmetric: PD.Boolean(false, { description: 'If true the automatic range is determined as [-|max|, |max|].' })
                })
            }),
            list: PD.ColorList('red-white-blue', { presetKind: 'scale' })
        })
    }),
    defaultColor: PD.Color(Color(0xcccccc)),
    normalOffset: PD.Numeric(0., { min: 0, max: 20, step: 0.1 }, { description: 'Offset vertex position along its normal by given amount.' }),
    usePalette: PD.Boolean(false, { description: 'Use a palette to color at the pixel level.' }),
};
export type ExternalVolumeColorThemeParams = typeof ExternalVolumeColorThemeParams

export function ExternalVolumeColorTheme(ctx: ThemeDataContext, props: PD.Values<ExternalVolumeColorThemeParams>): ColorTheme<ExternalVolumeColorThemeParams> {
    let volume: Volume | undefined;
    try {
        volume = props.volume.getValue();
    } catch {
        // .getValue() is resolved during state reconciliation => would throw from UI
    }

    // NOTE: this will currently be slow for with GPU/texture meshes due to slow iteration
    // TODO: create texture to be able to do the sampling on the GPU

    let color: LocationColor;
    let palette: ColorTheme.Palette | undefined;

    const { normalOffset, defaultColor, usePalette } = props;

    if (volume) {
        const coloring = props.coloring.params;
        const { stats } = volume.grid;
        const domain: [number, number] = coloring.domain.name === 'custom' ? coloring.domain.params : [stats.min, stats.max];

        const isRelative = props.coloring.name === 'relative-value';
        if (coloring.domain.name === 'auto' && isRelative) {
            domain[0] = (domain[0] - stats.mean) / stats.sigma;
            domain[1] = (domain[1] - stats.mean) / stats.sigma;
        }

        if (coloring.domain.name === 'auto' && coloring.domain.params.symmetric) {
            const max = Math.max(Math.abs(domain[0]), Math.abs(domain[1]));
            domain[0] = -max;
            domain[1] = max;
        }

        const scale = ColorScale.create({ domain, listOrName: coloring.list.colors });

        const position = Vec3();
        const getTrilinearlyInterpolated = Grid.makeGetTrilinearlyInterpolated(volume.grid, isRelative ? 'relative' : 'none');

        color = (location: Location): Color => {
            if (!isPositionLocation(location)) {
                return defaultColor;
            }

            // Offset the vertex position along its normal
            if (normalOffset > 0) {
                Vec3.scaleAndAdd(position, location.position, location.normal, normalOffset);
            } else {
                Vec3.copy(position, location.position);
            }

            const value = getTrilinearlyInterpolated(position);
            if (isNaN(value)) return defaultColor;

            if (usePalette) {
                return (clamp((value - domain[0]) / (domain[1] - domain[0]), 0, 1) * ColorTheme.PaletteScale) as Color;
            } else {
                return scale.color(value);
            }
        };

        palette = usePalette ? {
            colors: coloring.list.colors.map(e => Array.isArray(e) ? e[0] : e),
            filter: (coloring.list.kind === 'set' ? 'nearest' : 'linear') as 'nearest' | 'linear'
        } : undefined;
    } else {
        color = () => defaultColor;
    }

    return {
        factory: ExternalVolumeColorTheme,
        granularity: 'vertex',
        preferSmoothing: true,
        color,
        palette,
        props,
        description: Description,
        // TODO: figure out how to do legend for this
    };
}

export const ExternalVolumeColorThemeProvider: ColorTheme.Provider<ExternalVolumeColorThemeParams, 'external-volume'> = {
    name: 'external-volume',
    label: 'External Volume',
    category: ColorThemeCategory.Misc,
    factory: ExternalVolumeColorTheme,
    getParams: () => ExternalVolumeColorThemeParams,
    defaultValues: PD.getDefaultValues(ExternalVolumeColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => true,
};