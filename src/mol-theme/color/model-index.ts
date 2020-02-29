/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { StructureElement, Bond } from '../../mol-model/structure';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { ThemeDataContext } from '../../mol-theme/theme';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives every model a unique color based on the position (index) of the model in the list of models in the structure.'

export const ModelIndexColorThemeParams = {
    ...getPaletteParams({ type: 'scale', scaleList: 'red-yellow-blue' }),
}
export type ModelIndexColorThemeParams = typeof ModelIndexColorThemeParams
export function getModelIndexColorThemeParams(ctx: ThemeDataContext) {
    return ModelIndexColorThemeParams // TODO return copy
}

export function ModelIndexColorTheme(ctx: ThemeDataContext, props: PD.Values<ModelIndexColorThemeParams>): ColorTheme<ModelIndexColorThemeParams> {
    let color: LocationColor
    let legend: ScaleLegend | TableLegend | undefined

    if (ctx.structure) {
        const { models } = ctx.structure.root
        const palette = getPalette(models.length, props)
        legend = palette.legend
        const modelColor = new Map<string, Color>()
        for (let i = 0, il = models.length; i <il; ++i) {
            modelColor.set(models[i].id, palette.color(i))
        }

        color = (location: Location): Color => {
            if (StructureElement.Location.is(location)) {
                return modelColor.get(location.unit.model.id)!
            } else if (Bond.isLocation(location)) {
                return modelColor.get(location.aUnit.model.id)!
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        factory: ModelIndexColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    }
}

export const ModelIndexColorThemeProvider: ColorTheme.Provider<ModelIndexColorThemeParams> = {
    label: 'Model Index',
    category: ColorTheme.Category.Chain,
    factory: ModelIndexColorTheme,
    getParams: getModelIndexColorThemeParams,
    defaultValues: PD.getDefaultValues(ModelIndexColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.models.length > 1
}