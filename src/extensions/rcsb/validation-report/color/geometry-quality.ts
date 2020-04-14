/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ThemeDataContext } from '../../../../mol-theme/theme';
import { ColorTheme, LocationColor } from '../../../../mol-theme/color';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { Color } from '../../../../mol-util/color';
import { StructureElement } from '../../../../mol-model/structure';
import { Location } from '../../../../mol-model/location';
import { CustomProperty } from '../../../../mol-model-props/common/custom-property';
import { ValidationReportProvider, ValidationReport } from '../prop';
import { TableLegend } from '../../../../mol-util/legend';
import { PolymerType } from '../../../../mol-model/structure/model/types';
import { SetUtils } from '../../../../mol-util/set';

const DefaultColor = Color(0x909090);

const NoIssuesColor = Color(0x2166ac);
const OneIssueColor = Color(0xfee08b);
const TwoIssuesColor = Color(0xf46d43);
const ThreeOrMoreIssuesColor = Color(0xa50026);

const ColorLegend = TableLegend([
    ['Data unavailable', DefaultColor],
    ['No issues', NoIssuesColor],
    ['One issue', OneIssueColor],
    ['Two issues', TwoIssuesColor],
    ['Three or more issues', ThreeOrMoreIssuesColor],
]);

export function getGeometricQualityColorThemeParams(ctx: ThemeDataContext) {
    const validationReport = ctx.structure && ValidationReportProvider.get(ctx.structure.models[0]).value;
    const options: [string, string][] = [];
    if (validationReport) {
        const kinds = new Set<string>();
        validationReport.geometryIssues.forEach(v => v.forEach(k => kinds.add(k)));
        kinds.forEach(k => options.push([k, k]));
    }
    return {
        ignore: PD.MultiSelect([] as string[], options)
    };
}
export type GeometricQualityColorThemeParams = ReturnType<typeof getGeometricQualityColorThemeParams>

export function GeometryQualityColorTheme(ctx: ThemeDataContext, props: PD.Values<GeometricQualityColorThemeParams>): ColorTheme<GeometricQualityColorThemeParams> {
    let color: LocationColor = () => DefaultColor;

    const validationReport = ctx.structure && ValidationReportProvider.get(ctx.structure.models[0]);
    const contextHash = validationReport?.version;

    const value = validationReport?.value;
    const model = ctx.structure?.models[0];

    if (value && model) {
        const { geometryIssues, clashes, bondOutliers, angleOutliers } = value;
        const residueIndex = model.atomicHierarchy.residueAtomSegments.index;
        const { polymerType } = model.atomicHierarchy.derived.residue;
        const ignore = new Set(props.ignore);

        color = (location: Location): Color => {
            if (StructureElement.Location.is(location) && location.unit.model === model) {
                const { element } = location;
                const rI = residueIndex[element];

                const value = geometryIssues.get(rI);
                if (value === undefined) return DefaultColor;

                let count = SetUtils.differenceSize(value, ignore);

                if (count > 0 && polymerType[rI] === PolymerType.NA) {
                    count = 0;
                    if (!ignore.has('clash') && clashes.getVertexEdgeCount(element) > 0) count += 1;
                    if (!ignore.has('mog-bond-outlier') && bondOutliers.index.has(element)) count += 1;
                    if (!ignore.has('mog-angle-outlier') && angleOutliers.index.has(element)) count += 1;
                }

                switch (count) {
                    case undefined: return DefaultColor;
                    case 0: return NoIssuesColor;
                    case 1: return OneIssueColor;
                    case 2: return TwoIssuesColor;
                    default: return ThreeOrMoreIssuesColor;
                }
            }
            return DefaultColor;
        };
    }

    return {
        factory: GeometryQualityColorTheme,
        granularity: 'group',
        color,
        props,
        contextHash,
        description: 'Assigns residue colors according to the number of (filtered) geometry issues. Data from wwPDB Validation Report, obtained via RCSB PDB.',
        legend: ColorLegend
    };
}

export const GeometryQualityColorThemeProvider: ColorTheme.Provider<GeometricQualityColorThemeParams, ValidationReport.Tag.GeometryQuality> = {
    name: ValidationReport.Tag.GeometryQuality,
    label: 'Geometry Quality',
    category: ColorTheme.Category.Validation,
    factory: GeometryQualityColorTheme,
    getParams: getGeometricQualityColorThemeParams,
    defaultValues: PD.getDefaultValues(getGeometricQualityColorThemeParams({})),
    isApplicable: (ctx: ThemeDataContext) => ValidationReport.isApplicable(ctx.structure?.models[0]),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? ValidationReportProvider.attach(ctx, data.structure.models[0], void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && data.structure.models[0].customProperties.reference(ValidationReportProvider.descriptor, false)
    }
};