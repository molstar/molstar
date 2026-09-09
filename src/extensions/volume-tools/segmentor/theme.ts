/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ColorThemeCategory } from '../../../mol-theme/color/categories';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Color } from '../../../mol-util/color';
import { TableLegend } from '../../../mol-util/legend';
import { Location } from '../../../mol-model/location';
import { Grid, Volume } from '../../../mol-model/volume';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { isPositionLocation } from '../../../mol-geo/util/location-iterator';
import { BodyLabels } from './labels';
import { MaxBodyId } from './types';

const Description = 'Colors a volume surface by the body each voxel is assigned to.';

export const BodyLabelColorThemeParams = {
    unassignedColor: PD.Color(Color(0x9a9a9a), { description: 'Color of voxels not assigned to any body.' }),
    /** Mirrors `LabelStore.version`; bumping it forces a color update after labels change. */
    version: PD.Numeric(0, {}, { isHidden: true }),
};
export type BodyLabelColorThemeParams = typeof BodyLabelColorThemeParams

function clampCell(v: number, n: number) {
    return Math.min(Math.max(v, 0), Math.max(n - 2, 0));
}

/**
 * Maps a world position to the label of the densest corner of the grid cell containing it.
 * Isosurface vertices lie on cell edges between an inside and an outside voxel; taking the
 * densest corner always picks the inside voxel, so surfaces are colored without speckle.
 */
export function makeLabelAtPosition(volume: Volume, labels: Uint8Array): (position: Vec3) => number {
    const { space, data } = volume.grid.cells;
    const values = data as unknown as ArrayLike<number>;
    const [nx, ny, nz] = space.dimensions as [number, number, number];
    const c2g = Mat4.invert(Mat4(), Grid.getGridToCartesianTransform(volume.grid));
    const g = Vec3();

    return (position: Vec3) => {
        Vec3.transformMat4(g, position, c2g);
        const i0 = clampCell(Math.floor(g[0]), nx), i1 = Math.min(i0 + 1, nx - 1);
        const j0 = clampCell(Math.floor(g[1]), ny), j1 = Math.min(j0 + 1, ny - 1);
        const k0 = clampCell(Math.floor(g[2]), nz), k1 = Math.min(k0 + 1, nz - 1);

        let best = -Infinity, label = 0;
        for (let i = i0; i <= i1; i++) {
            for (let j = j0; j <= j1; j++) {
                for (let k = k0; k <= k1; k++) {
                    const o = space.dataOffset(i, j, k);
                    const v = values[o];
                    if (v > best) {
                        best = v;
                        label = labels[o];
                    }
                }
            }
        }
        return label;
    };
}

export function BodyLabelColorTheme(ctx: ThemeDataContext, props: PD.Values<BodyLabelColorThemeParams>): ColorTheme<BodyLabelColorThemeParams> {
    const volume = ctx.volume;
    const store = volume && BodyLabels.get(volume);

    if (!volume || !store) {
        return {
            factory: BodyLabelColorTheme,
            granularity: 'uniform',
            color: () => props.unassignedColor,
            props,
            description: Description,
        };
    }

    const colors = new Array<Color>(MaxBodyId + 1).fill(props.unassignedColor);
    for (const body of store.bodies) colors[body.id] = body.color;
    const labelAt = makeLabelAtPosition(volume, store.labels);

    const color: LocationColor = (location: Location) => {
        if (!isPositionLocation(location)) return props.unassignedColor;
        return colors[labelAt(location.position)];
    };

    return {
        factory: BodyLabelColorTheme,
        granularity: 'vertex',
        preferSmoothing: false,
        color,
        props,
        description: Description,
        contextHash: store.version,
        legend: TableLegend(store.bodies.map(b => [b.name, b.color] as [string, Color])),
    };
}

export const BodyLabelColorThemeProvider: ColorTheme.Provider<BodyLabelColorThemeParams, 'body-label'> = {
    name: 'body-label',
    label: 'Body Label',
    category: ColorThemeCategory.Misc,
    factory: BodyLabelColorTheme,
    getParams: () => BodyLabelColorThemeParams,
    defaultValues: PD.getDefaultValues(BodyLabelColorThemeParams),
    // Applicable to any volume so the theme survives param normalisation before labels exist.
    isApplicable: (ctx: ThemeDataContext) => !!ctx.volume,
};
