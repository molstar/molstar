/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { ErrorColor, NtCColors } from '@molstar/dnatco-extension/color';
import { ConfalPyramidsProvider } from './property.js';
import { ConfalPyramidsTypes as CPT } from './types.js';
import { Dnatco } from '@molstar/dnatco-extension/property';
import type { Location } from '@molstar/model/model/location';
import { CustomProperty } from '@molstar/model/props/common/custom-property';
import { ColorTheme } from '@molstar/graphics/theme/color';
import type { ThemeDataContext } from '@molstar/graphics/theme/theme';
import { Color, ColorMap } from '@molstar/core/util/color';
import { getColorMapParams } from '@molstar/core/util/color/params';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { TableLegend } from '@molstar/core/util/legend';
import { ObjectKeys } from '@molstar/core/util/type-helpers';
import { ColorThemeCategory } from '@molstar/graphics/theme/color/categories';

const Description = 'Assigns colors to confal pyramids';

const PyramidsColors = ColorMap({ ...NtCColors });
type PyramidsColors = typeof PyramidsColors;

export const ConfalPyramidsColorThemeParams = {
    colors: PD.MappedStatic('default', {
        'default': PD.EmptyGroup(),
        'custom': PD.Group(getColorMapParams(PyramidsColors))
    }),
};
export type ConfalPyramidsColorThemeParams = typeof ConfalPyramidsColorThemeParams;

export function getConfalPyramidsColorThemeParams(ctx: ThemeDataContext) {
    return PD.clone(ConfalPyramidsColorThemeParams);
}

export function ConfalPyramidsColorTheme(ctx: ThemeDataContext, props: PD.Values<ConfalPyramidsColorThemeParams>): ColorTheme<ConfalPyramidsColorThemeParams> {
    const colorMap = props.colors.name === 'default' ? PyramidsColors : props.colors.params;

    function color(location: Location, isSecondary: boolean): Color {
        if (CPT.isLocation(location)) {
            const { step, isLower } = location.data;
            const key = step.NtC + `_${isLower ? 'Lwr' : 'Upr'}` as keyof PyramidsColors;
            return colorMap[key] ?? ErrorColor;
        }

        return ErrorColor;
    }

    return {
        factory: ConfalPyramidsColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: TableLegend(ObjectKeys(colorMap).map(k => [k.replace('_', ' '), colorMap[k]] as [string, Color]).concat([['Error', ErrorColor]])),
    };
}

export const ConfalPyramidsColorThemeProvider: ColorTheme.Provider<ConfalPyramidsColorThemeParams, 'confal-pyramids'> = {
    name: 'confal-pyramids',
    label: 'Confal Pyramids',
    category: ColorThemeCategory.Residue,
    factory: ConfalPyramidsColorTheme,
    getParams: getConfalPyramidsColorThemeParams,
    defaultValues: PD.getDefaultValues(ConfalPyramidsColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.models.some(m => Dnatco.isApplicable(m)),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? ConfalPyramidsProvider.attach(ctx, data.structure.models[0], void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && ConfalPyramidsProvider.ref(data.structure.models[0], false)
    }
};
