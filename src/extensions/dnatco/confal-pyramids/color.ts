/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { ConfalPyramidsProvider } from './property';
import { ConfalPyramidsTypes as CPT } from './types';
import { DnatcoCommon as DC } from '../common';
import { Location } from '../../../mol-model/location';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { ColorTheme } from '../../../mol-theme/color';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Assigns colors to confal pyramids';
const ErrorColor = Color(0xFFA10A);

function getConformerColor(ntc: string, useLower: boolean, palette: PD.Values<ConfalPyramidsColorThemeParams>): Color {
    const item = CPT.NtCToClasses.get(ntc);
    if (!item) return ErrorColor;

    const key = useLower ? item[1] : item[0];
    if (!palette.hasOwnProperty(key)) return ErrorColor;
    return palette[key];
}

export const ConfalPyramidsColorThemeParams = {
    A: PD.Color(Color(0xFFC1C1)),
    B: PD.Color(Color(0xC8CFFF)),
    BII: PD.Color(Color(0x0059DA)),
    miB: PD.Color(Color(0x3BE8FB)),
    Z: PD.Color(Color(0x01F60E)),
    IC: PD.Color(Color(0xFA5CFB)),
    OPN: PD.Color(Color(0xE90000)),
    SYN: PD.Color(Color(0xFFFF01)),
    N: PD.Color(Color(0xF2F2F2))
};
export type ConfalPyramidsColorThemeParams = typeof ConfalPyramidsColorThemeParams
export function getConfalPyramidsColorThemeParams(ctx: ThemeDataContext) {
    return PD.clone(ConfalPyramidsColorThemeParams);
}

export function ConfalPyramidsColorTheme(ctx: ThemeDataContext, props: PD.Values<ConfalPyramidsColorThemeParams>): ColorTheme<ConfalPyramidsColorThemeParams> {
    function color(location: Location, isSecondary: boolean): Color {
        if (CPT.isLocation(location)) {
            const { pyramid, isLower } = location.data;
            return getConformerColor(pyramid.NtC, isLower, props);
        }

        return DefaultColor;
    }

    return {
        factory: ConfalPyramidsColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
    };
}

export const ConfalPyramidsColorThemeProvider: ColorTheme.Provider<ConfalPyramidsColorThemeParams, 'confal-pyramids'> = {
    name: 'confal-pyramids',
    label: 'Confal Pyramids',
    category: ColorTheme.Category.Residue,
    factory: ConfalPyramidsColorTheme,
    getParams: getConfalPyramidsColorThemeParams,
    defaultValues: PD.getDefaultValues(ConfalPyramidsColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.models.some(m => DC.isApplicable(m)),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? ConfalPyramidsProvider.attach(ctx, data.structure.models[0], void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && data.structure.models[0].customProperties.reference(ConfalPyramidsProvider.descriptor, false)
    }
};
