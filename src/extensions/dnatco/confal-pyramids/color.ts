/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { ConfalPyramids, ConfalPyramidsProvider } from './property';
import { ConfalPyramidsTypes as CPT } from './types';
import { Location } from '../../../mol-model/location';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { ColorTheme } from '../../../mol-theme/color';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { Color, ColorMap } from '../../../mol-util/color';
import { getColorMapParams } from '../../../mol-util/color/params';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { TableLegend } from '../../../mol-util/legend';
import { ObjectKeys } from '../../../mol-util/type-helpers';

const Description = 'Assigns colors to confal pyramids';

const DefaultClassColors = {
    A: 0xFFC1C1,
    B: 0xC8CFFF,
    BII: 0x0059DA,
    miB: 0x3BE8FB,
    Z: 0x01F60E,
    IC: 0xFA5CFB,
    OPN: 0xE90000,
    SYN: 0xFFFF01,
    N: 0xF2F2F2,
};
const ErrorColor = Color(0xFFA10A);

const PyramidsColors = ColorMap({
    NANT_Upr: DefaultClassColors.N,
    NANT_Lwr: DefaultClassColors.N,
    AA00_Upr: DefaultClassColors.A,
    AA00_Lwr: DefaultClassColors.A,
    AA02_Upr: DefaultClassColors.A,
    AA02_Lwr: DefaultClassColors.A,
    AA03_Upr: DefaultClassColors.A,
    AA03_Lwr: DefaultClassColors.A,
    AA04_Upr: DefaultClassColors.A,
    AA04_Lwr: DefaultClassColors.A,
    AA08_Upr: DefaultClassColors.A,
    AA08_Lwr: DefaultClassColors.A,
    AA09_Upr: DefaultClassColors.A,
    AA09_Lwr: DefaultClassColors.A,
    AA01_Upr: DefaultClassColors.A,
    AA01_Lwr: DefaultClassColors.A,
    AA05_Upr: DefaultClassColors.A,
    AA05_Lwr: DefaultClassColors.A,
    AA06_Upr: DefaultClassColors.A,
    AA06_Lwr: DefaultClassColors.A,
    AA10_Upr: DefaultClassColors.A,
    AA10_Lwr: DefaultClassColors.A,
    AA11_Upr: DefaultClassColors.A,
    AA11_Lwr: DefaultClassColors.A,
    AA07_Upr: DefaultClassColors.A,
    AA07_Lwr: DefaultClassColors.A,
    AA12_Upr: DefaultClassColors.A,
    AA12_Lwr: DefaultClassColors.A,
    AA13_Upr: DefaultClassColors.A,
    AA13_Lwr: DefaultClassColors.A,
    AB01_Upr: DefaultClassColors.A,
    AB01_Lwr: DefaultClassColors.B,
    AB02_Upr: DefaultClassColors.A,
    AB02_Lwr: DefaultClassColors.B,
    AB03_Upr: DefaultClassColors.A,
    AB03_Lwr: DefaultClassColors.B,
    AB04_Upr: DefaultClassColors.A,
    AB04_Lwr: DefaultClassColors.B,
    AB05_Upr: DefaultClassColors.A,
    AB05_Lwr: DefaultClassColors.B,
    BA01_Upr: DefaultClassColors.B,
    BA01_Lwr: DefaultClassColors.A,
    BA05_Upr: DefaultClassColors.B,
    BA05_Lwr: DefaultClassColors.A,
    BA09_Upr: DefaultClassColors.B,
    BA09_Lwr: DefaultClassColors.A,
    BA08_Upr: DefaultClassColors.BII,
    BA08_Lwr: DefaultClassColors.A,
    BA10_Upr: DefaultClassColors.B,
    BA10_Lwr: DefaultClassColors.A,
    BA13_Upr: DefaultClassColors.BII,
    BA13_Lwr: DefaultClassColors.A,
    BA16_Upr: DefaultClassColors.BII,
    BA16_Lwr: DefaultClassColors.A,
    BA17_Upr: DefaultClassColors.BII,
    BA17_Lwr: DefaultClassColors.A,
    BB00_Upr: DefaultClassColors.B,
    BB00_Lwr: DefaultClassColors.B,
    BB01_Upr: DefaultClassColors.B,
    BB01_Lwr: DefaultClassColors.B,
    BB17_Upr: DefaultClassColors.B,
    BB17_Lwr: DefaultClassColors.B,
    BB02_Upr: DefaultClassColors.B,
    BB02_Lwr: DefaultClassColors.B,
    BB03_Upr: DefaultClassColors.B,
    BB03_Lwr: DefaultClassColors.B,
    BB11_Upr: DefaultClassColors.B,
    BB11_Lwr: DefaultClassColors.B,
    BB16_Upr: DefaultClassColors.B,
    BB16_Lwr: DefaultClassColors.B,
    BB04_Upr: DefaultClassColors.B,
    BB04_Lwr: DefaultClassColors.BII,
    BB05_Upr: DefaultClassColors.B,
    BB05_Lwr: DefaultClassColors.BII,
    BB07_Upr: DefaultClassColors.BII,
    BB07_Lwr: DefaultClassColors.BII,
    BB08_Upr: DefaultClassColors.BII,
    BB08_Lwr: DefaultClassColors.BII,
    BB10_Upr: DefaultClassColors.miB,
    BB10_Lwr: DefaultClassColors.miB,
    BB12_Upr: DefaultClassColors.miB,
    BB12_Lwr: DefaultClassColors.miB,
    BB13_Upr: DefaultClassColors.miB,
    BB13_Lwr: DefaultClassColors.miB,
    BB14_Upr: DefaultClassColors.miB,
    BB14_Lwr: DefaultClassColors.miB,
    BB15_Upr: DefaultClassColors.miB,
    BB15_Lwr: DefaultClassColors.miB,
    BB20_Upr: DefaultClassColors.miB,
    BB20_Lwr: DefaultClassColors.miB,
    IC01_Upr: DefaultClassColors.IC,
    IC01_Lwr: DefaultClassColors.IC,
    IC02_Upr: DefaultClassColors.IC,
    IC02_Lwr: DefaultClassColors.IC,
    IC03_Upr: DefaultClassColors.IC,
    IC03_Lwr: DefaultClassColors.IC,
    IC04_Upr: DefaultClassColors.IC,
    IC04_Lwr: DefaultClassColors.IC,
    IC05_Upr: DefaultClassColors.IC,
    IC05_Lwr: DefaultClassColors.IC,
    IC06_Upr: DefaultClassColors.IC,
    IC06_Lwr: DefaultClassColors.IC,
    IC07_Upr: DefaultClassColors.IC,
    IC07_Lwr: DefaultClassColors.IC,
    OP01_Upr: DefaultClassColors.OPN,
    OP01_Lwr: DefaultClassColors.OPN,
    OP02_Upr: DefaultClassColors.OPN,
    OP02_Lwr: DefaultClassColors.OPN,
    OP03_Upr: DefaultClassColors.OPN,
    OP03_Lwr: DefaultClassColors.OPN,
    OP04_Upr: DefaultClassColors.OPN,
    OP04_Lwr: DefaultClassColors.OPN,
    OP05_Upr: DefaultClassColors.OPN,
    OP05_Lwr: DefaultClassColors.OPN,
    OP06_Upr: DefaultClassColors.OPN,
    OP06_Lwr: DefaultClassColors.OPN,
    OP07_Upr: DefaultClassColors.OPN,
    OP07_Lwr: DefaultClassColors.OPN,
    OP08_Upr: DefaultClassColors.OPN,
    OP08_Lwr: DefaultClassColors.OPN,
    OP09_Upr: DefaultClassColors.OPN,
    OP09_Lwr: DefaultClassColors.OPN,
    OP10_Upr: DefaultClassColors.OPN,
    OP10_Lwr: DefaultClassColors.OPN,
    OP11_Upr: DefaultClassColors.OPN,
    OP11_Lwr: DefaultClassColors.OPN,
    OP12_Upr: DefaultClassColors.OPN,
    OP12_Lwr: DefaultClassColors.OPN,
    OP13_Upr: DefaultClassColors.OPN,
    OP13_Lwr: DefaultClassColors.OPN,
    OP14_Upr: DefaultClassColors.OPN,
    OP14_Lwr: DefaultClassColors.OPN,
    OP15_Upr: DefaultClassColors.OPN,
    OP15_Lwr: DefaultClassColors.OPN,
    OP16_Upr: DefaultClassColors.OPN,
    OP16_Lwr: DefaultClassColors.OPN,
    OP17_Upr: DefaultClassColors.OPN,
    OP17_Lwr: DefaultClassColors.OPN,
    OP18_Upr: DefaultClassColors.OPN,
    OP18_Lwr: DefaultClassColors.OPN,
    OP19_Upr: DefaultClassColors.OPN,
    OP19_Lwr: DefaultClassColors.OPN,
    OP20_Upr: DefaultClassColors.OPN,
    OP20_Lwr: DefaultClassColors.OPN,
    OP21_Upr: DefaultClassColors.OPN,
    OP21_Lwr: DefaultClassColors.OPN,
    OP22_Upr: DefaultClassColors.OPN,
    OP22_Lwr: DefaultClassColors.OPN,
    OP23_Upr: DefaultClassColors.OPN,
    OP23_Lwr: DefaultClassColors.OPN,
    OP24_Upr: DefaultClassColors.OPN,
    OP24_Lwr: DefaultClassColors.OPN,
    OP25_Upr: DefaultClassColors.OPN,
    OP25_Lwr: DefaultClassColors.OPN,
    OP26_Upr: DefaultClassColors.OPN,
    OP26_Lwr: DefaultClassColors.OPN,
    OP27_Upr: DefaultClassColors.OPN,
    OP27_Lwr: DefaultClassColors.OPN,
    OP28_Upr: DefaultClassColors.OPN,
    OP28_Lwr: DefaultClassColors.OPN,
    OP29_Upr: DefaultClassColors.OPN,
    OP29_Lwr: DefaultClassColors.OPN,
    OP30_Upr: DefaultClassColors.OPN,
    OP30_Lwr: DefaultClassColors.OPN,
    OP31_Upr: DefaultClassColors.OPN,
    OP31_Lwr: DefaultClassColors.OPN,
    OPS1_Upr: DefaultClassColors.OPN,
    OPS1_Lwr: DefaultClassColors.OPN,
    OP1S_Upr: DefaultClassColors.OPN,
    OP1S_Lwr: DefaultClassColors.SYN,
    AAS1_Upr: DefaultClassColors.SYN,
    AAS1_Lwr: DefaultClassColors.A,
    AB1S_Upr: DefaultClassColors.A,
    AB1S_Lwr: DefaultClassColors.SYN,
    AB2S_Upr: DefaultClassColors.A,
    AB2S_Lwr: DefaultClassColors.SYN,
    BB1S_Upr: DefaultClassColors.B,
    BB1S_Lwr: DefaultClassColors.SYN,
    BB2S_Upr: DefaultClassColors.B,
    BB2S_Lwr: DefaultClassColors.SYN,
    BBS1_Upr: DefaultClassColors.SYN,
    BBS1_Lwr: DefaultClassColors.B,
    ZZ01_Upr: DefaultClassColors.Z,
    ZZ01_Lwr: DefaultClassColors.Z,
    ZZ02_Upr: DefaultClassColors.Z,
    ZZ02_Lwr: DefaultClassColors.Z,
    ZZ1S_Upr: DefaultClassColors.Z,
    ZZ1S_Lwr: DefaultClassColors.SYN,
    ZZ2S_Upr: DefaultClassColors.Z,
    ZZ2S_Lwr: DefaultClassColors.SYN,
    ZZS1_Upr: DefaultClassColors.SYN,
    ZZS1_Lwr: DefaultClassColors.Z,
    ZZS2_Upr: DefaultClassColors.SYN,
    ZZS2_Lwr: DefaultClassColors.Z,
});
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
            const { pyramid, isLower } = location.data;
            const key = pyramid.NtC + `_${isLower ? 'Lwr' : 'Upr'}` as keyof PyramidsColors;
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
    category: ColorTheme.Category.Residue,
    factory: ConfalPyramidsColorTheme,
    getParams: getConfalPyramidsColorThemeParams,
    defaultValues: PD.getDefaultValues(ConfalPyramidsColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.models.some(m => ConfalPyramids.isApplicable(m)),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => data.structure ? ConfalPyramidsProvider.attach(ctx, data.structure.models[0], void 0, true) : Promise.resolve(),
        detach: (data) => data.structure && ConfalPyramidsProvider.ref(data.structure.models[0], false)
    }
};
