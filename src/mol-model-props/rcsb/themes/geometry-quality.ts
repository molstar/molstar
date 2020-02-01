/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ThemeDataContext } from '../../../mol-theme/theme';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition'
import { Color } from '../../../mol-util/color';
import { StructureElement } from '../../../mol-model/structure';
import { Location } from '../../../mol-model/location';
import { CustomProperty } from '../../common/custom-property';
import { ValidationReportProvider, ValidationReport } from '../validation-report';
import { TableLegend } from '../../../mol-util/legend';
import { PolymerType } from '../../../mol-model/structure/model/types';

const DefaultColor = Color(0x909090)

const NoIssuesColor = Color(0x2166ac)
const OneIssueColor = Color(0xfee08b)
const TwoIssuesColor = Color(0xf46d43)
const ThreeOrMoreIssuesColor = Color(0xa50026)

const ColorLegend = TableLegend([
    ['No issues', NoIssuesColor],
    ['One issue', OneIssueColor],
    ['Two issues', TwoIssuesColor],
    ['Three or more issues', ThreeOrMoreIssuesColor]
])

export function GeometryQualityColorTheme(ctx: ThemeDataContext, props: {}): ColorTheme<{}> {
    let color: LocationColor = () => DefaultColor

    const validationReport = ctx.structure && ValidationReportProvider.get(ctx.structure.models[0])
    const contextHash = validationReport?.version

    const value = validationReport?.value
    const model = ctx.structure?.models[0]

    if (value && model) {
        const { geometryIssues, clashes, bondOutliers, angleOutliers } = value
        const residueIndex = model.atomicHierarchy.residueAtomSegments.index
        const { polymerType } = model.atomicHierarchy.derived.residue
        color = (location: Location): Color => {
            if (StructureElement.Location.is(location) && location.unit.model === model) {
                const { element } = location
                const rI = residueIndex[element]
                let value = geometryIssues.get(rI)?.size
                if (value !== undefined && polymerType[rI] === PolymerType.NA) {
                    value = 0
                    if (clashes.getVertexEdgeCount(element) > 0) value += 1
                    if (bondOutliers.index.has(element)) value += 1
                    if (angleOutliers.index.has(element)) value += 1
                }

                switch (value) {
                    case undefined: return DefaultColor
                    case 0: return NoIssuesColor
                    case 1: return OneIssueColor
                    case 2: return TwoIssuesColor
                    default: return ThreeOrMoreIssuesColor
                }
            }
            return DefaultColor
        }
    }

    return {
        factory: GeometryQualityColorTheme,
        granularity: 'group',
        color,
        props,
        contextHash,
        description: 'Assigns residue colors according to the number of geometry issues.',
        legend: ColorLegend
    }
}

export const GeometryQualityColorThemeProvider: ColorTheme.Provider<{}> = {
    label: 'RCSB Geometry Quality',
    factory: GeometryQualityColorTheme,
    getParams: () => ({}),
    defaultValues: PD.getDefaultValues({}),
    isApplicable: (ctx: ThemeDataContext) => ValidationReport.isApplicable(ctx.structure?.models[0]),
    ensureCustomProperties: (ctx: CustomProperty.Context, data: ThemeDataContext) => {
        return data.structure ? ValidationReportProvider.attach(ctx, data.structure.models[0]) : Promise.resolve()
    }
}