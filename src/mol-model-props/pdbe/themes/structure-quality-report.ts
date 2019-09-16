/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureQualityReport } from '../../../mol-model-props/pdbe/structure-quality-report';
import { Location } from '../../../mol-model/location';
import { StructureElement } from '../../../mol-model/structure';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { Color } from '../../../mol-util/color';
import { TableLegend } from '../../../mol-util/legend';

const ValidationColors = [
    Color.fromRgb(170, 170, 170), // not applicable
    Color.fromRgb(0, 255, 0), // 0 issues
    Color.fromRgb(255, 255, 0), // 1
    Color.fromRgb(255, 128, 0), // 2
    Color.fromRgb(255, 0, 0), // 3 or more
]

const ValidationColorTable: [string, Color][] = [
    ['No Issues', ValidationColors[1]],
    ['One Issue', ValidationColors[2]],
    ['Two Issues', ValidationColors[3]],
    ['Three Or More Issues', ValidationColors[4]],
    ['Not Applicable', ValidationColors[9]]
]

export function StructureQualityReportColorTheme(ctx: ThemeDataContext, props: {}): ColorTheme<{}> {
    let color: LocationColor

    if (ctx.structure && !ctx.structure.isEmpty && ctx.structure.models[0].customProperties.has(StructureQualityReport.Descriptor)) {
        const getIssues = StructureQualityReport.getIssues;
        color = (location: Location) => {
            if (StructureElement.Location.is(location)) {
                return ValidationColors[Math.min(3, getIssues(location).length) + 1];
            }
            return ValidationColors[0];
        }
    } else {
        color = () => ValidationColors[0];
    }

    return {
        factory: StructureQualityReportColorTheme,
        granularity: 'group',
        color: color,
        props: props,
        description: 'Assigns residue colors according to the number of issues in the PDBe Validation Report.',
        legend: TableLegend(ValidationColorTable)
    }
}