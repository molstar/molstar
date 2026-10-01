/**
 * Copyright (c) 2021-25 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { QualityAssessment, QualityAssessmentProvider } from '../prop.js';
import type { Location } from '@molstar/model/model/location';
import { Bond, StructureElement, Unit } from '@molstar/model/model/structure';
import { ColorTheme, type LocationColor } from '@molstar/graphics/theme/color';
import type { ThemeDataContext } from '@molstar/graphics/theme/theme';
import { Color, ColorScale } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperty } from '@molstar/model/props/common/custom-property';
import { ColorThemeCategory } from '@molstar/graphics/theme/color/categories';

const DefaultColor = Color(0xaaaaaa);

export function getQmeanScoreColorThemeParams(ctx: ThemeDataContext) {
    return {
        metricId: QualityAssessment.getLocalOptions(ctx.structure?.models[0], 'qmean'),
    };
}
export type QmeanScoreColorThemeParams = ReturnType<typeof getQmeanScoreColorThemeParams>

export function QmeanScoreColorTheme(ctx: ThemeDataContext, props: PD.Values<QmeanScoreColorThemeParams>): ColorTheme<QmeanScoreColorThemeParams> {
    let color: LocationColor = () => DefaultColor;

    const scale = ColorScale.create({
        domain: [0, 1],
        listOrName: [
            [Color(0xFF5000), 0.5], [Color(0x025AFD), 1.0]
        ]
    });

    if (ctx.structure) {
        const l = StructureElement.Location.create(ctx.structure.root);

        const getColor = (location: StructureElement.Location): Color => {
            const { unit, element } = location;
            if (!Unit.isAtomic(unit)) return DefaultColor;
            const qualityAssessment = QualityAssessmentProvider.get(unit.model).value;
            const metric = qualityAssessment?.localMap.get(props.metricId!)?.values ?? qualityAssessment?.qmean;
            const score = metric?.get(unit.model.atomicHierarchy.residueAtomSegments.index[element]) ?? -1;
            if (score < 0) {
                return DefaultColor;
            } else {
                return scale.color(score);
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
        factory: QmeanScoreColorTheme,
        granularity: 'group',
        preferSmoothing: true,
        color,
        props,
        description: 'Assigns residue colors according to the QMEAN score.',
        legend: scale.legend
    };
}

export const QmeanScoreColorThemeProvider: ColorTheme.Provider<QmeanScoreColorThemeParams, 'qmean-score'> = {
    name: 'qmean-score',
    label: 'QMEAN Score',
    category: ColorThemeCategory.Validation,
    factory: QmeanScoreColorTheme,
    getParams: getQmeanScoreColorThemeParams,
    defaultValues: PD.getDefaultValues(getQmeanScoreColorThemeParams({})),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure?.models.some(m => QualityAssessment.isApplicable(m, 'qmean')),
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