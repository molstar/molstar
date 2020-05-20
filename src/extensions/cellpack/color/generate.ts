/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ThemeDataContext } from '../../../mol-theme/theme';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Color } from '../../../mol-util/color';
import { getPalette } from '../../../mol-util/color/palette';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ScaleLegend, TableLegend } from '../../../mol-util/legend';
import { StructureElement, Bond, Model } from '../../../mol-model/structure';
import { Location } from '../../../mol-model/location';
import { CellPackInfoProvider } from '../property';
import { distinctColors } from '../../../mol-util/color/distinct';
import { Hcl } from '../../../mol-util/color/spaces/hcl';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every model in a CellPack packing a unique generated color similar to other models in the packing.';

export const CellPackGenerateColorThemeParams = {};
export type CellPackGenerateColorThemeParams = typeof CellPackGenerateColorThemeParams
export function getCellPackGenerateColorThemeParams(ctx: ThemeDataContext) {
    return CellPackGenerateColorThemeParams; // TODO return copy
}

export function CellPackGenerateColorTheme(ctx: ThemeDataContext, props: PD.Values<CellPackGenerateColorThemeParams>): ColorTheme<CellPackGenerateColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    const info = ctx.structure && CellPackInfoProvider.get(ctx.structure).value;

    if (ctx.structure && info) {
        const colors = distinctColors(info.packingsCount);
        let hcl = Hcl.fromColor(Hcl(), colors[info.packingIndex]);

        const hue = [Math.max(0, hcl[0] - 35), Math.min(360, hcl[0] + 35)] as [number, number];

        const { models } = ctx.structure.root;

        let size = 0;
        for (const m of models) size = Math.max(size, Model.TrajectoryInfo.get(m).size);

        const palette = getPalette(size, { palette: {
            name: 'generate',
            params: {
                hue, chroma: [30, 80], luminance: [15, 85],
                clusteringStepCount: 50, minSampleCount: 800, maxCount: 75,
                minLabel: 'Min', maxLabel: 'Max', valueLabel: (i: number) => `${i + 1}`,
            }
        }});
        legend = palette.legend;
        const modelColor = new Map<number, Color>();
        for (let i = 0, il = models.length; i < il; ++i) {
            const idx = Model.TrajectoryInfo.get(models[i]).index;
            modelColor.set(Model.TrajectoryInfo.get(models[i]).index, palette.color(idx));
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
        factory: CellPackGenerateColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    };
}

export const CellPackGenerateColorThemeProvider: ColorTheme.Provider<CellPackGenerateColorThemeParams, 'cellpack-generate'> = {
    name: 'cellpack-generate',
    label: 'CellPack Generate',
    category: ColorTheme.Category.Chain,
    factory: CellPackGenerateColorTheme,
    getParams: getCellPackGenerateColorThemeParams,
    defaultValues: PD.getDefaultValues(CellPackGenerateColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => {
        return (
            !!ctx.structure && ctx.structure.elementCount > 0 &&
            Model.TrajectoryInfo.get(ctx.structure.models[0]).size > 1 &&
            !!CellPackInfoProvider.get(ctx.structure).value
        );
    }
};

