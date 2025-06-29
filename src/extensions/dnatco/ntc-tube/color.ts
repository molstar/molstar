/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { ErrorColor, NtCColors } from '../color';
import { NtCTubeProvider } from './property';
import { NtCTubeTypes as NTT } from './types';
import { Dnatco } from '../property';
import { Location } from '../../../mol-model/location';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { ColorTheme } from '../../../mol-theme/color';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { Color, ColorMap } from '../../../mol-util/color';
import { getColorMapParams } from '../../../mol-util/color/params';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { TableLegend } from '../../../mol-util/legend';
import { ObjectKeys } from '../../../mol-util/type-helpers';
import { ColorThemeCategory } from '../../../mol-theme/color/categories';

const Description = 'Assigns colors to NtC Tube segments';

const NtCTubeColors = ColorMap({
    ...NtCColors,
    residueMarker: Color(0x222222),
    stepBoundaryMarker: Color(0x656565),
});
type NtCTubeColors = typeof NtCTubeColors;

export const NtCTubeColorThemeParams = {
    colors: PD.MappedStatic('default', {
        'default': PD.EmptyGroup(),
        'custom': PD.Group(getColorMapParams(NtCTubeColors)),
        'uniform': PD.Color(Color(0xEEEEEE)),
    }),
    markResidueBoundaries: PD.Boolean(true),
    markSegmentBoundaries: PD.Boolean(true),
};
export type NtCTubeColorThemeParams = typeof NtCTubeColorThemeParams;

export function getNtCTubeColorThemeParams(ctx: ThemeDataContext) {
    return PD.clone(NtCTubeColorThemeParams);
}

export function NtCTubeColorTheme(ctx: ThemeDataContext, props: PD.Values<NtCTubeColorThemeParams>): ColorTheme<NtCTubeColorThemeParams> {
    const colorMap = props.colors.name === 'default'
        ? NtCTubeColors
        : props.colors.name === 'custom'
            ? props.colors.params
            : ColorMap({
                ...Object.fromEntries(ObjectKeys(NtCTubeColors).map(item => [item, props.colors.params])),
                residueMarker: NtCTubeColors.residueMarker,
                stepBoundaryMarker: NtCTubeColors.stepBoundaryMarker
            }) as NtCTubeColors;

    function color(location: Location, isSecondary: boolean): Color {
        if (NTT.isLocation(location)) {
            const { data } = location;
            const { step, kind } = data;
            let key;
            if (kind === 'upper')
                key = step.NtC + '_Upr' as keyof NtCTubeColors;
            else if (kind === 'lower')
                key = step.NtC + '_Lwr' as keyof NtCTubeColors;
            else if (kind === 'residue-boundary')
                key = (!props.markResidueBoundaries ? step.NtC + '_Lwr' : 'residueMarker') as keyof NtCTubeColors;
            else /* segment-boundary */
                key = (!props.markSegmentBoundaries ? step.NtC + '_Lwr' : 'stepBoundaryMarker') as keyof NtCTubeColors;

            return colorMap[key] ?? ErrorColor;
        }

        return ErrorColor;
    }

    return {
        factory: NtCTubeColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: TableLegend(ObjectKeys(colorMap).map(k => [k.replace('_', ' '), colorMap[k]] as [string, Color]).concat([['Error', ErrorColor]])),
    };
}

export const NtCTubeColorThemeProvider: ColorTheme.Provider<NtCTubeColorThemeParams, 'ntc-tube'> = {
    name: 'ntc-tube',
    label: 'NtC Tube',
    category: ColorThemeCategory.Residue,
    factory: NtCTubeColorTheme,
    getParams: getNtCTubeColorThemeParams,
    defaultValues: PD.getDefaultValues(NtCTubeColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.models.some(m => Dnatco.isApplicable(m)),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? NtCTubeProvider.attach(ctx, data.structure.models[0], void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && NtCTubeProvider.ref(data.structure.models[0], false)
    }
};
