/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

const Description = 'Assign size based on the given value of a volume cell. Negative values are made positive.';

export const VolumeValueSizeThemeParams = {
    scale: PD.Numeric(1, { min: 0.1, max: 5, step: 0.1 }),
};
export type VolumeValueSizeThemeParams = typeof VolumeValueSizeThemeParams
export function getVolumeValueSizeThemeParams(ctx: ThemeDataContext) {
    return VolumeValueSizeThemeParams; // TODO return copy
}

export function VolumeValueSizeTheme(ctx: ThemeDataContext, props: PD.Values<VolumeValueSizeThemeParams>): SizeTheme<VolumeValueSizeThemeParams> {
    if (ctx.volume) {
        const { data } = ctx.volume.grid.cells;

        const isLocation = Volume.Cell.isLocation;
        const size: LocationSize = (location: Location): number => {
            if (isLocation(location)) {
                return Math.abs(data[location.cell]) * props.scale;
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
