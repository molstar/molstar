/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { StructureElement, Bond, Model } from '../../mol-model/structure';
import type { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';
import { ColorThemeCategory } from './categories';
import { hashFnv32a } from '../../mol-data/util';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every model (frame) a unique color based on the index in its trajectory.';

export const TrajectoryIndexColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: 'purples' }),
};
export type TrajectoryIndexColorThemeParams = typeof TrajectoryIndexColorThemeParams
export function getTrajectoryIndexColorThemeParams(ctx: ThemeDataContext) {
    return PD.clone(TrajectoryIndexColorThemeParams);
}

export function TrajectoryIndexColorTheme(ctx: ThemeDataContext, props: PD.Values<TrajectoryIndexColorThemeParams>): ColorTheme<TrajectoryIndexColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;
    let contextHash = -1;

    if (ctx.structure) {
        const { models } = ctx.structure.root;

        let size = 0;
        for (const m of models) size = Math.max(size, Model.TrajectoryInfo.get(m)?.size || 0);

        const hash: number[] = [size];
        const palette = getPalette(size, props);
        legend = palette.legend;
        const modelColor = new Map<number, Color>();
        for (let i = 0, il = models.length; i < il; ++i) {
            const idx = Model.TrajectoryInfo.get(models[i])?.index || 0;
            modelColor.set(idx, palette.color(idx));
            hash.push(idx);
        }
        contextHash = hashFnv32a(hash);

        color = (location: Location): Color => {
            if (StructureElement.Location.is(location)) {
                return modelColor.get(Model.TrajectoryInfo.get(location.unit.model).index)!;
            } else if (Bond.isLocation(location)) {
                return modelColor.get(Model.TrajectoryInfo.get(location.aUnit.model).index)!;
            }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: TrajectoryIndexColorTheme,
        granularity: 'instance',
        color,
        props,
        contextHash,
        description: Description,
        legend
    };
}

export const TrajectoryIndexColorThemeProvider: ColorTheme.Provider<TrajectoryIndexColorThemeParams, 'trajectory-index'> = {
    name: 'trajectory-index',
    label: 'Trajectory Index',
    category: ColorThemeCategory.Chain,
    factory: TrajectoryIndexColorTheme,
    getParams: getTrajectoryIndexColorThemeParams,
    defaultValues: PD.getDefaultValues(TrajectoryIndexColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.elementCount > 0 && Model.TrajectoryInfo.get(ctx.structure.models[0]).size > 1
};