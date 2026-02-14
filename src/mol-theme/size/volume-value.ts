/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { SizeTheme } from '../size';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { LocationSize } from '../../mol-geo/geometry/size-data';
import { Volume } from '../../mol-model/volume/volume';
import { Location } from '../../mol-model/location';
import { isPositionLocation } from '../../mol-geo/util/location-iterator';
import { Grid } from '../../mol-model/volume/grid';
import { clamp } from '../../mol-math/interpolate';

const Description = 'Assign size based on the given value of a volume cell.';

export const VolumeValueSizeThemeParams = {
    scale: PD.Numeric(1, { min: 0.01, max: 10, step: 0.01 }),
    transform: PD.Select('linear', [['linear', 'Linear'], ['quadratic', 'Quadratic'], ['cubic', 'Cubic']] as const, { description: 'Linear: value * scale, Quadratic: value * scale^2, Cubic: value * scale^3' }),
    domain: PD.MappedStatic('auto', {
        custom: PD.Interval([0, 1], { step: 0.001, min: 0 }),
        auto: PD.Group({})
    }),
};
export type VolumeValueSizeThemeParams = typeof VolumeValueSizeThemeParams
export function getVolumeValueSizeThemeParams(ctx: ThemeDataContext) {
    return VolumeValueSizeThemeParams; // TODO return copy
}

export function VolumeValueSizeTheme(ctx: ThemeDataContext, props: PD.Values<VolumeValueSizeThemeParams>): SizeTheme<VolumeValueSizeThemeParams> {
    if (ctx.volume) {
        const { min, max } = ctx.volume.grid.stats;
        const domain: [number, number] = props.domain.name === 'custom' ? props.domain.params : [min, max];
        const scaleFactor = props.transform === 'cubic' ? props.scale ** 3 : props.transform === 'quadratic' ? props.scale ** 2 : props.scale;

        if (ctx.locationKinds?.includes('cell-location')) {
            const { data } = ctx.volume.grid.cells;

            const isLocation = Volume.Cell.isLocation;
            const size: LocationSize = (location: Location): number => {
                if (isLocation(location)) {
                    const v = clamp(Math.abs(data[location.cell]), domain[0], domain[1]);
                    return v * scaleFactor;
                } else {
                    return 0;
                }
            };

            return {
                factory: VolumeValueSizeTheme,
                granularity: 'group',
                size,
                props,
                description: Description
            };
        } else {
            const getTrilinearlyInterpolated = Grid.makeGetTrilinearlyInterpolated(ctx.volume.grid, 'none');

            const size: LocationSize = (location: Location): number => {
                if (isPositionLocation(location)) {
                    const value = getTrilinearlyInterpolated(location.position);
                    if (Number.isNaN(value)) return 0;
                    const v = clamp(Math.abs(value), domain[0], domain[1]);
                    return v * scaleFactor;
                } else {
                    return 0;
                }
            };

            return {
                factory: VolumeValueSizeTheme,
                granularity: 'vertex',
                size,
                props,
                description: Description
            };
        }
    } else {
        return {
            factory: VolumeValueSizeTheme,
            granularity: 'uniform',
            size: () => props.scale,
            props,
            description: Description
        };
    }
}

export const VolumeValueSizeThemeProvider: SizeTheme.Provider<VolumeValueSizeThemeParams, 'volume-value'> = {
    name: 'volume-value',
    label: 'Volume Value',
    category: '',
    factory: VolumeValueSizeTheme,
    getParams: getVolumeValueSizeThemeParams,
    defaultValues: PD.getDefaultValues(VolumeValueSizeThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.volume && !Volume.Segmentation.get(ctx.volume),
};
