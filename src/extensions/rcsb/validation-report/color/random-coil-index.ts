/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ThemeDataContext } from '../../../../mol-theme/theme';
import { ColorTheme, LocationColor } from '../../../../mol-theme/color';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { Color, ColorScale } from '../../../../mol-util/color';
import { StructureElement, Model, ElementIndex, Bond } from '../../../../mol-model/structure';
import { Location } from '../../../../mol-model/location';
import { CustomProperty } from '../../../../mol-model-props/common/custom-property';
import { ValidationReportProvider, ValidationReport } from '../prop';
import { ColorThemeCategory } from '../../../../mol-theme/color/categories';

const DefaultColor = Color(0xCCCCCC);

export function RandomCoilIndexColorTheme(ctx: ThemeDataContext, props: {}): ColorTheme<{}> {
    let color: LocationColor = () => DefaultColor;
    const scale = ColorScale.create({
        reverse: true,
        domain: [0, 0.6],
        listOrName: 'red-yellow-blue',
    });

    const validationReport = ctx.structure && ValidationReportProvider.get(ctx.structure.models[0]);
    const contextHash = validationReport?.version;

    const rci = validationReport?.value?.rci;
    const model = ctx.structure?.models[0];

    if (rci && model) {
        const residueIndex = model.atomicHierarchy.residueAtomSegments.index;
        const getColor = (element: ElementIndex) => {
            const value = rci.get(residueIndex[element]);
            return value === undefined ? DefaultColor : scale.color(value);
        };

        color = (location: Location): Color => {
            if (StructureElement.Location.is(location) && location.unit.model === model) {
                return getColor(location.element);
            } else if (Bond.isLocation(location) && location.aUnit.model === model) {
                return getColor(location.aUnit.elements[location.aIndex]);
            }
            return DefaultColor;
        };
    }

    return {
        factory: RandomCoilIndexColorTheme,
        granularity: 'group',
        preferSmoothing: true,
        color,
        props,
        contextHash,
        description: 'Assigns residue colors according to the Random Coil Index value. Data from wwPDB Validation Report, obtained via RCSB PDB.',
        legend: scale.legend
    };
}

export const RandomCoilIndexColorThemeProvider: ColorTheme.Provider<{}, ValidationReport.Tag.RandomCoilIndex> = {
    name: ValidationReport.Tag.RandomCoilIndex,
    label: 'Random Coil Index',
    category: ColorThemeCategory.Validation,
    factory: RandomCoilIndexColorTheme,
    getParams: () => ({}),
    defaultValues: PD.getDefaultValues({}),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ValidationReport.isApplicable(ctx.structure.models[0]) && Model.isFromNmr(ctx.structure.models[0]),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? ValidationReportProvider.attach(ctx, data.structure.models[0], void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && ValidationReportProvider.ref(data.structure.models[0], false)
    }
};