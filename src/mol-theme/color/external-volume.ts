/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Cai Huiyu <szmun.caihy@gmail.com>
 */

import { Color, ColorScale } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import type { ColorTheme } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { Grid, Volume } from '../../mol-model/volume';
import { type PluginContext } from '../../mol-plugin/context';
import { isPositionLocation } from '../../mol-geo/util/location-iterator';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { lerp } from '../../mol-math/interpolate';
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
                custom: PD.Interval([-1, 1]),
                auto: PD.Group({
                    symmetric: PD.Boolean(false, { description: 'If true the automatic range is determined as [-|max|, |max|].' })
                })
            }),
            list: PD.ColorList('red-white-blue', { presetKind: 'scale' })
        }),
        'relative-value': PD.Group({
            domain: PD.MappedStatic('auto', {
                custom: PD.Interval([-1, 1]),
                auto: PD.Group({
                    symmetric: PD.Boolean(false, { description: 'If true the automatic range is determined as [-|max|, |max|].' })
                })
            }),
            list: PD.ColorList('red-white-blue', { presetKind: 'scale' })
        })
    }),
    defaultColor: PD.Color(Color(0xcccccc)),
    normalOffset: PD.Numeric(0., { min: 0, max: 20, step: 0.1 }, { description: 'Offset vertex position along its normal by given amount.' }),
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

    let color;
    if (volume) {
        const coloring = props.coloring.params;
        const { stats } = volume.grid;
        const domain: [number, number] = coloring.domain.name === 'custom' ? coloring.domain.params : [stats.min, stats.max];

        const isRelative = props.coloring.name === 'relative-value';
        if (coloring.domain.name === 'auto' && isRelative) {
            domain[0] = (domain[0] - stats.mean) / stats.sigma;
            domain[1] = (domain[1] - stats.mean) / stats.sigma;
        }

        if (props.coloring.params.domain.name === 'auto' && props.coloring.params.domain.params.symmetric) {
            const max = Math.max(Math.abs(domain[0]), Math.abs(domain[1]));
            domain[0] = -max;
            domain[1] = max;
        }

        const scale = ColorScale.create({ domain, listOrName: coloring.list.colors });

        const cartnToGrid = Grid.getGridToCartesianTransform(volume.grid);
        Mat4.invert(cartnToGrid, cartnToGrid);
        const gridCoords = Vec3();

        const { dimensions, get } = volume.grid.cells.space;
        const data = volume.grid.cells.data;

        const [mi, mj, mk] = dimensions;

        color = (location: Location): Color => {
            if (!isPositionLocation(location)) {
                return props.defaultColor;
            }

            // Offset the vertex position along its normal
            Vec3.copy(gridCoords, location.position);

            if (props.normalOffset > 0) {
                Vec3.scaleAndAdd(gridCoords, gridCoords, location.normal, props.normalOffset);
            }

            Vec3.transformMat4(gridCoords, gridCoords, cartnToGrid);

            const i = Math.floor(gridCoords[0]);
            const j = Math.floor(gridCoords[1]);
            const k = Math.floor(gridCoords[2]);

            if (i < 0 || i >= mi || j < 0 || j >= mj || k < 0 || k >= mk) {
                return props.defaultColor;
            }

            const u = gridCoords[0] - i;
            const v = gridCoords[1] - j;
            const w = gridCoords[2] - k;

            // Tri-linear interpolation for the value
            const ii = Math.min(i + 1, mi - 1);
            const jj = Math.min(j + 1, mj - 1);
            const kk = Math.min(k + 1, mk - 1);

            let a = get(data, i, j, k);
            let b = get(data, ii, j, k);
            let c = get(data, i, jj, k);
            let d = get(data, ii, jj, k);
            const x = lerp(lerp(a, b, u), lerp(c, d, u), v);

            a = get(data, i, j, kk);
            b = get(data, ii, j, kk);
            c = get(data, i, jj, kk);
            d = get(data, ii, jj, kk);
            const y = lerp(lerp(a, b, u), lerp(c, d, u), v);

            let value = lerp(x, y, w);
            if (isRelative) {
                value = (value - stats.mean) / stats.sigma;
            }

            return scale.color(value);
        };
    } else {
        color = () => props.defaultColor;
    }

    return {
        factory: ExternalVolumeColorTheme,
        granularity: 'vertex',
        preferSmoothing: true,
        color,
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