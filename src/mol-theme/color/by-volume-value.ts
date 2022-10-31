/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color, ColorScale } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { ColorTheme } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { Grid, Volume } from '../../mol-model/volume';
import { type PluginContext } from '../../mol-plugin/context';
import { isPositionLocation } from '../../mol-geo/util/location-iterator';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { lerp } from '../../mol-math/interpolate';

const Description = `Assigns a color based volume value at a given vertex.`;

export const ByVolumeValueColorThemeParams = {
    volume: PD.ValueRef<Volume>(
        (ctx: PluginContext) => {
            const volumes = ctx.state.data.selectQ(q => q.root.subtree().filter(c => Volume.is(c.obj?.data)));
            return volumes.map(v => [v.transform.ref, v.obj?.label ?? '<unknown>'] as [string, string]);
        },
        (ref, getData) => getData(ref),
    ),
    domain: PD.MappedStatic('auto', {
        custom: PD.Interval([-1, 1]),
        auto: PD.EmptyGroup()
    }),
    list: PD.ColorList('red-white-blue', { presetKind: 'scale' }),
    defaultColor: PD.Color(Color(0xcccccc))
};
export type ByVolumeValueColorThemeParams = typeof ByVolumeValueColorThemeParams

export function ByVolumeValueColorTheme(ctx: ThemeDataContext, props: PD.Values<ByVolumeValueColorThemeParams>): ColorTheme<ByVolumeValueColorThemeParams> {
    let volume: Volume | undefined;
    try {
        volume = props.volume.getValue();
    } catch {
        // .getValue() is resolved during state reconciliation => would throw from UI
    }

    const domain: [number, number] = props.domain.name === 'custom' ? props.domain.params : volume ? [volume.grid.stats.min, volume.grid.stats.max] : [-1, 1];
    const scale = ColorScale.create({ domain, listOrName: props.list.colors });

    // NOTE: this will currently not work with GPU iso-surfaces since it requires vertex coloring
    // TODO: create texture to be able to do the sampling on the GPU

    let color;
    if (volume) {
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

            Vec3.copy(gridCoords, location.position);
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

            const value = lerp(x, y, w);

            return scale.color(value);
        };
    } else {
        color = () => props.defaultColor;
    }

    return {
        factory: ByVolumeValueColorTheme,
        granularity: 'vertex',
        preferSmoothing: true,
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    };
}

export const ByVolumeValueColorThemeProvider: ColorTheme.Provider<ByVolumeValueColorThemeParams, 'by-volume-value'> = {
    name: 'by-volume-value',
    label: 'By Volume Value',
    category: ColorTheme.Category.Misc,
    factory: ByVolumeValueColorTheme,
    getParams: () => ByVolumeValueColorThemeParams,
    defaultValues: PD.getDefaultValues(ByVolumeValueColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => true, // TODO
};