/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { StructureElement, Bond, Model } from '../../mol-model/structure';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every model a unique color based on the position (index) of the model in the list of models in the structure.';

export const ModelIndexColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: 'purples' }),
};
export type ModelIndexColorThemeParams = typeof ModelIndexColorThemeParams
export function getModelIndexColorThemeParams(ctx: ThemeDataContext) {
    return ModelIndexColorThemeParams; // TODO return copy
}

export function ModelIndexColorTheme(ctx: ThemeDataContext, props: PD.Values<ModelIndexColorThemeParams>): ColorTheme<ModelIndexColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    if (ctx.structure) {
        const { models } = ctx.structure.root;

        let size = 0;
        for (const m of models) size = Math.max(size, Model.TrajectoryInfo.get(m)?.size || 0);

        const palette = getPalette(size, props);
        legend = palette.legend;
        const modelColor = new Map<number, Color>();
        for (let i = 0, il = models.length; i < il; ++i) {
            const idx = Model.TrajectoryInfo.get(models[i])?.index || 0;
            modelColor.set(idx, palette.color(idx));
        }

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
        factory: ModelIndexColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    };
}

export const ModelIndexColorThemeProvider: ColorTheme.Provider<ModelIndexColorThemeParams, 'model-index'> = {
    name: 'model-index',
    label: 'Model Index',
    category: ColorTheme.Category.Chain,
    factory: ModelIndexColorTheme,
    getParams: getModelIndexColorThemeParams,
    defaultValues: PD.getDefaultValues(ModelIndexColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.elementCount > 0 && Model.TrajectoryInfo.get(ctx.structure.models[0]).size > 1
};