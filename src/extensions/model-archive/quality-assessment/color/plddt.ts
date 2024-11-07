/**
 * Copyright (c) 2021-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Mandar Deshpande <mandar@ebi.ac.uk>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QualityAssessment, QualityAssessmentProvider } from '../prop';
import { Location } from '../../../../mol-model/location';
import { Bond, Model, StructureElement, Unit } from '../../../../mol-model/structure';
import { ColorTheme, LocationColor } from '../../../../mol-theme/color';
import { ThemeDataContext } from '../../../../mol-theme/theme';
import { Color } from '../../../../mol-util/color';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { CustomProperty } from '../../../../mol-model-props/common/custom-property';
import { TableLegend } from '../../../../mol-util/legend';

const DefaultColor = Color(0xaaaaaa);
const ConfidenceColors = {
    'No Score': DefaultColor,
    'Very Low': Color(0xff7d45),
    'Low': Color(0xffdb13),
    'Confident': Color(0x65cbf3),
    'Very High': Color(0x0053d6)
};

const ConfidenceColorLegend = TableLegend(Object.entries(ConfidenceColors));

export function getPLDDTConfidenceColorThemeParams(ctx: ThemeDataContext) {
    return {};
}
export type PLDDTConfidenceColorThemeParams = ReturnType<typeof getPLDDTConfidenceColorThemeParams>

export function PLDDTConfidenceColorTheme(ctx: ThemeDataContext, props: PD.Values<PLDDTConfidenceColorThemeParams>): ColorTheme<PLDDTConfidenceColorThemeParams> {
    let color: LocationColor = () => DefaultColor;

    if (ctx.structure) {
        const l = StructureElement.Location.create(ctx.structure.root);

        const getColor = (location: StructureElement.Location): Color => {
            const { unit, element } = location;
            if (!Unit.isAtomic(unit)) return DefaultColor;

            const qualityAssessment = QualityAssessmentProvider.get(unit.model).value;
            let score = qualityAssessment?.pLDDT?.get(unit.model.atomicHierarchy.residueAtomSegments.index[element]);
            if (typeof score !== 'number') {
                score = unit.model.atomicConformation.B_iso_or_equiv.value(element);
            }

            if (score < 0) {
                return DefaultColor;
            } else if (score <= 50) {
                return Color(0xff7d45);
            } else if (score <= 70) {
                return Color(0xffdb13);
            } else if (score <= 90) {
                return Color(0x65cbf3);
            } else {
                return Color(0x0053d6);
            }
        };

        color = (location: Location) => {
            if (StructureElement.Location.is(location)) {
                return getColor(location);
            } else if (Bond.isLocation(location)) {
                l.unit = location.aUnit;
                l.element = location.aUnit.elements[location.aIndex];
                return getColor(l);
            }
            return DefaultColor;
        };
    }

    return {
        factory: PLDDTConfidenceColorTheme,
        granularity: 'group',
        preferSmoothing: true,
        color,
        props,
        description: 'Assigns residue colors according to the pLDDT Confidence score. If no Model Archive quality assessment score is available, the B-factor value is used instead.',
        legend: ConfidenceColorLegend
    };
}

export const PLDDTConfidenceColorThemeProvider: ColorTheme.Provider<PLDDTConfidenceColorThemeParams, 'plddt-confidence'> = {
    name: 'plddt-confidence',
    label: 'pLDDT Confidence',
    category: ColorTheme.Category.Validation,
    factory: PLDDTConfidenceColorTheme,
    getParams: getPLDDTConfidenceColorThemeParams,
    defaultValues: PD.getDefaultValues(getPLDDTConfidenceColorThemeParams({})),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure?.models.some(m => QualityAssessment.isApplicable(m, 'pLDDT') || (m.atomicConformation.B_iso_or_equiv.isDefined && !Model.isExperimental(m))),
    ensureCustomProperties: {
        attach: async (ctx: CustomProperty.Context, data: ThemeDataContext) => {
            if (data.structure) {
                for (const m of data.structure.models) {
                    await QualityAssessmentProvider.attach(ctx, m, void 0, true);
                }
            }
        },
        detach: async (data: ThemeDataContext) => {
            if (data.structure) {
                for (const m of data.structure.models) {
                    QualityAssessmentProvider.ref(m, false);
                }
            }
        }
    }
};