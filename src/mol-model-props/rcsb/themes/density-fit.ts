/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ThemeDataContext } from '../../../mol-theme/theme';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition'
import { Color, ColorScale } from '../../../mol-util/color';
import { StructureElement } from '../../../mol-model/structure';
import { Location } from '../../../mol-model/location';
import { CustomProperty } from '../../common/custom-property';
import { ValidationReportProvider, ValidationReport } from '../validation-report';

const DefaultColor = Color(0xCCCCCC)

export function DensityFitColorTheme(ctx: ThemeDataContext, props: {}): ColorTheme<{}> {
    let color: LocationColor = () => DefaultColor
    const scaleRsrz = ColorScale.create({
        minLabel: 'Poor',
        maxLabel: 'Better',
        domain: [2, 0],
        listOrName: 'red-yellow-blue',
    })
    const scaleRscc = ColorScale.create({
        minLabel: 'Poor',
        maxLabel: 'Better',
        domain: [0.678, 1.0],
        listOrName: 'red-yellow-blue',
    })

    const validationReport = ctx.structure && ValidationReportProvider.get(ctx.structure.models[0])
    const contextHash = validationReport?.version
    const model = ctx.structure?.models[0]

    if (validationReport?.value && model) {
        const { rsrz, rscc } = validationReport.value
        const residueIndex = model.atomicHierarchy.residueAtomSegments.index
        color = (location: Location): Color => {
            if (StructureElement.Location.is(location) && location.unit.model === model) {
                const rsrzValue = rsrz.get(residueIndex[location.element])
                if (rsrzValue !== undefined) return scaleRsrz.color(rsrzValue)
                const rsccValue = rscc.get(residueIndex[location.element])
                if (rsccValue !== undefined) return scaleRscc.color(rsccValue)
                return DefaultColor
            }
            return DefaultColor
        }
    }

    return {
        factory: DensityFitColorTheme,
        granularity: 'group',
        color,
        props,
        contextHash,
        description: 'Assigns residue colors according to the density fit using normalized Real Space R (RSRZ) for polymer residues and real space correlation coefficient (RSCC) for ligands. Colors range from poor (RSRZ = 2 or RSCC = 0.678) - to better (RSRZ = 0 or RSCC = 1.0). Data from wwPDB Validation Report, obtained via RCSB PDB.',
        legend: scaleRsrz.legend
    }
}

export const DensityFitColorThemeProvider: ColorTheme.Provider<{}> = {
    label: 'Density Fit',
    category: ColorTheme.Category.Validation,
    factory: DensityFitColorTheme,
    getParams: () => ({}),
    defaultValues: PD.getDefaultValues({}),
    isApplicable: (ctx: ThemeDataContext) => ValidationReport.isApplicable(ctx.structure?.models[0]),
    ensureCustomProperties: (ctx: CustomProperty.Context, data: ThemeDataContext) => {
        return data.structure ? ValidationReportProvider.attach(ctx, data.structure.models[0]) : Promise.resolve()
    }
}