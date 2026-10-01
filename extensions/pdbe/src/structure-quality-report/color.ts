/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureQualityReport, StructureQualityReportProvider } from './prop.js';
import type { Location } from '@molstar/model/model/location';
import { Bond, StructureElement } from '@molstar/model/model/structure';
import { ColorTheme, type LocationColor } from '@molstar/graphics/theme/color';
import type { ThemeDataContext } from '@molstar/graphics/theme/theme';
import { Color } from '@molstar/core/util/color';
import { TableLegend } from '@molstar/core/util/legend';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperty } from '@molstar/model/props/common/custom-property';
import { ColorThemeCategory } from '@molstar/graphics/theme/color/categories';

const ValidationColors = [
    Color.fromRgb(170, 170, 170), // not applicable
    Color.fromRgb(0, 255, 0), // 0 issues
    Color.fromRgb(255, 255, 0), // 1
    Color.fromRgb(255, 128, 0), // 2
    Color.fromRgb(255, 0, 0), // 3 or more
];

const ValidationColorTable: [string, Color][] = [
    ['No Issues', ValidationColors[1]],
    ['One Issue', ValidationColors[2]],
    ['Two Issues', ValidationColors[3]],
    ['Three Or More Issues', ValidationColors[4]],
    ['Not Applicable', ValidationColors[0]],
];

export const StructureQualityReportColorThemeParams = {
    type: PD.MappedStatic('issue-count', {
        'issue-count': PD.Group({}),
        'specific-issue': PD.Group({
            kind: PD.Text()
        })
    })
};

type Params = typeof StructureQualityReportColorThemeParams

export function StructureQualityReportColorTheme(ctx: ThemeDataContext, props: PD.Values<Params>): ColorTheme<Params> {
    let color: LocationColor;

    if (ctx.structure && !ctx.structure.isEmpty && ctx.structure.models[0].customProperties.has(StructureQualityReportProvider.descriptor)) {
        const getIssues = StructureQualityReport.getIssues;
        const l = StructureElement.Location.create(ctx.structure);

        if (props.type.name === 'issue-count') {
            color = (location: Location) => {
                if (StructureElement.Location.is(location)) {
                    return ValidationColors[Math.min(3, getIssues(location).length) + 1];
                } else if (Bond.isLocation(location)) {
                    l.unit = location.aUnit;
                    l.element = location.aUnit.elements[location.aIndex];
                    return ValidationColors[Math.min(3, getIssues(l).length) + 1];
                }
                return ValidationColors[0];
            };
        } else {
            const issue = props.type.params.kind;
            color = (location: Location) => {
                if (StructureElement.Location.is(location) && getIssues(location).indexOf(issue) >= 0) {
                    return ValidationColors[4];
                } else if (Bond.isLocation(location)) {
                    l.unit = location.aUnit;
                    l.element = location.aUnit.elements[location.aIndex];
                    return ValidationColors[Math.min(3, getIssues(l).length) + 1];
                }
                return ValidationColors[0];
            };
        }
    } else {
        color = () => ValidationColors[0];
    }

    return {
        factory: StructureQualityReportColorTheme,
        granularity: 'group',
        preferSmoothing: true,
        color: color,
        props: props,
        description: 'Assigns residue colors according to the number of quality issues or a specific quality issue. Data from wwPDB Validation Report, obtained via PDBe.',
        legend: TableLegend(ValidationColorTable)
    };
}

export const StructureQualityReportColorThemeProvider: ColorTheme.Provider<Params, 'pdbe-structure-quality-report'> = {
    name: 'pdbe-structure-quality-report',
    label: 'Structure Quality Report',
    category: ColorThemeCategory.Validation,
    factory: StructureQualityReportColorTheme,
    getParams: ctx => {
        const issueTypes = StructureQualityReport.getIssueTypes(ctx.structure);
        if (issueTypes.length === 0) {
            return {
                type: PD.MappedStatic('issue-count', {
                    'issue-count': PD.Group({})
                })
            };
        }

        return {
            type: PD.MappedStatic('issue-count', {
                'issue-count': PD.Group({}),
                'specific-issue': PD.Group({
                    kind: PD.Select(issueTypes[0], PD.arrayToOptions(issueTypes))
                }, { isFlat: true })
            })
        };
    },
    defaultValues: PD.getDefaultValues(StructureQualityReportColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => StructureQualityReport.isApplicable(ctx.structure?.models[0]),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? StructureQualityReportProvider.attach(ctx, data.structure.models[0], void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && StructureQualityReportProvider.ref(data.structure.models[0], false)
    }
};