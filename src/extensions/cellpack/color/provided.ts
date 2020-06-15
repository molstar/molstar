/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ThemeDataContext } from '../../../mol-theme/theme';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Color } from '../../../mol-util/color';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ScaleLegend, TableLegend } from '../../../mol-util/legend';
import { StructureElement, Model } from '../../../mol-model/structure';
import { Location } from '../../../mol-model/location';
import { CellPackInfoProvider } from '../property';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every model in a CellPack the color provied in the packing data.';

export const CellPackProvidedColorThemeParams = {};
export type CellPackProvidedColorThemeParams = typeof CellPackProvidedColorThemeParams
export function getCellPackProvidedColorThemeParams(ctx: ThemeDataContext) {
    return CellPackProvidedColorThemeParams; // TODO return copy
}

export function CellPackProvidedColorTheme(ctx: ThemeDataContext, props: PD.Values<CellPackProvidedColorThemeParams>): ColorTheme<CellPackProvidedColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    const info = ctx.structure && CellPackInfoProvider.get(ctx.structure).value;

    if (ctx.structure && info?.colors) {
        const { models } = ctx.structure.root;
        const modelColor = new Map<number, Color>();
        for (let i = 0, il = models.length; i < il; ++i) {
            const idx = Model.TrajectoryInfo.get(models[i]).index;
            modelColor.set(Model.TrajectoryInfo.get(models[i]).index, info.colors[idx]);
        }

        color = (location: Location): Color => {
            return StructureElement.Location.is(location)
                ? modelColor.get(Model.TrajectoryInfo.get(location.unit.model).index)!
                : DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: CellPackProvidedColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    };
}

export const CellPackProvidedColorThemeProvider: ColorTheme.Provider<CellPackProvidedColorThemeParams, 'cellpack-provided'> = {
    name: 'cellpack-provided',
    label: 'CellPack Provided',
    category: ColorTheme.Category.Chain,
    factory: CellPackProvidedColorTheme,
    getParams: getCellPackProvidedColorThemeParams,
    defaultValues: PD.getDefaultValues(CellPackProvidedColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => {
        return (
            !!ctx.structure && ctx.structure.elementCount > 0 &&
            Model.TrajectoryInfo.get(ctx.structure.models[0]).size > 1 &&
            !!CellPackInfoProvider.get(ctx.structure).value?.colors
        );
    }
};