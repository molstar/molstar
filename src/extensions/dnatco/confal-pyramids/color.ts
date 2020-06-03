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
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { TableLegend } from '../../../mol-util/legend';
import { iterableToArray } from '../../../mol-data/util';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Assigns colors to confal pyramids';
const ErrorColor = Color(0xFFA10A);

type ConformerClasses = 'A' | 'B' | 'BII' | 'miB' | 'Z' | 'IC' | 'OPN' | 'SYN' | 'N';

const ColorMapping: ReadonlyMap<ConformerClasses, Color> = new Map([
    ['A', Color(0xFFC1C1)],
    ['B', Color(0xC8CFFF)],
    ['BII', Color(0x0059DA)],
    ['miB', Color(0x3BE8FB)],
    ['Z',  Color(0x01F60E)],
    ['IC', Color(0xFA5CFB)],
    ['OPN', Color(0xE90000)],
    ['SYN', Color(0xFFFF01)],
    ['N', Color(0xF2F2F2)],
]);

const NtCToClasses: ReadonlyMap<string, [ConformerClasses, ConformerClasses]> = new Map([
    ['NANT', ['N', 'N']],
    ['AA00', ['A', 'A']],
    ['AA02', ['A', 'A']],
    ['AA03', ['A', 'A']],
    ['AA04', ['A', 'A']],
    ['AA08', ['A', 'A']],
    ['AA09', ['A', 'A']],
    ['AA01', ['A', 'A']],
    ['AA05', ['A', 'A']],
    ['AA06', ['A', 'A']],
    ['AA10', ['A', 'A']],
    ['AA11', ['A', 'A']],
    ['AA07', ['A', 'A']],
    ['AA12', ['A', 'A']],
    ['AA13', ['A', 'A']],
    ['AB01', ['A', 'B']],
    ['AB02', ['A', 'B']],
    ['AB03', ['A', 'B']],
    ['AB04', ['A', 'B']],
    ['AB05', ['A', 'B']],
    ['BA01', ['B', 'A']],
    ['BA05', ['B', 'A']],
    ['BA09', ['B', 'A']],
    ['BA08', ['BII', 'A']],
    ['BA10', ['B', 'A']],
    ['BA13', ['BII', 'A']],
    ['BA16', ['BII', 'A']],
    ['BA17', ['BII', 'A']],
    ['BB00', ['B', 'B']],
    ['BB01', ['B', 'B']],
    ['BB17', ['B', 'B']],
    ['BB02', ['B', 'B']],
    ['BB03', ['B', 'B']],
    ['BB11', ['B', 'B']],
    ['BB16', ['B', 'B']],
    ['BB04', ['B', 'BII']],
    ['BB05', ['B', 'BII']],
    ['BB07', ['BII', 'BII']],
    ['BB08', ['BII', 'BII']],
    ['BB10', ['miB', 'miB']],
    ['BB12', ['miB', 'miB']],
    ['BB13', ['miB', 'miB']],
    ['BB14', ['miB', 'miB']],
    ['BB15', ['miB', 'miB']],
    ['BB20', ['miB', 'miB']],
    ['IC01', ['IC', 'IC']],
    ['IC02', ['IC', 'IC']],
    ['IC03', ['IC', 'IC']],
    ['IC04', ['IC', 'IC']],
    ['IC05', ['IC', 'IC']],
    ['IC06', ['IC', 'IC']],
    ['IC07', ['IC', 'IC']],
    ['OP01', ['OPN', 'OPN']],
    ['OP02', ['OPN', 'OPN']],
    ['OP03', ['OPN', 'OPN']],
    ['OP04', ['OPN', 'OPN']],
    ['OP05', ['OPN', 'OPN']],
    ['OP06', ['OPN', 'OPN']],
    ['OP07', ['OPN', 'OPN']],
    ['OP08', ['OPN', 'OPN']],
    ['OP09', ['OPN', 'OPN']],
    ['OP10', ['OPN', 'OPN']],
    ['OP11', ['OPN', 'OPN']],
    ['OP12', ['OPN', 'OPN']],
    ['OP13', ['OPN', 'OPN']],
    ['OP14', ['OPN', 'OPN']],
    ['OP15', ['OPN', 'OPN']],
    ['OP16', ['OPN', 'OPN']],
    ['OP17', ['OPN', 'OPN']],
    ['OP18', ['OPN', 'OPN']],
    ['OP19', ['OPN', 'OPN']],
    ['OP20', ['OPN', 'OPN']],
    ['OP21', ['OPN', 'OPN']],
    ['OP22', ['OPN', 'OPN']],
    ['OP23', ['OPN', 'OPN']],
    ['OP24', ['OPN', 'OPN']],
    ['OP25', ['OPN', 'OPN']],
    ['OP26', ['OPN', 'OPN']],
    ['OP27', ['OPN', 'OPN']],
    ['OP28', ['OPN', 'OPN']],
    ['OP29', ['OPN', 'OPN']],
    ['OP30', ['OPN', 'OPN']],
    ['OP31', ['OPN', 'OPN']],
    ['OPS1', ['OPN', 'OPN']],
    ['OP1S', ['OPN', 'SYN']],
    ['AAS1', ['SYN', 'A']],
    ['AB1S', ['A', 'SYN']],
    ['AB2S', ['A', 'SYN']],
    ['BB1S', ['B', 'SYN']],
    ['BB2S', ['B', 'SYN']],
    ['BBS1', ['SYN', 'B']],
    ['ZZ01', ['Z', 'Z']],
    ['ZZ02', ['Z', 'Z']],
    ['ZZ1S', ['Z', 'SYN']],
    ['ZZ2S', ['Z', 'SYN']],
    ['ZZS1', ['SYN', 'Z']],
    ['ZZS2', ['SYN', 'Z']],
]);

function getConformerColor(ntc: string, useLower: boolean): Color {
    const item = NtCToClasses.get(ntc);
    if (!item) return ErrorColor;
    return ColorMapping.get(useLower ? item[1] : item[0]) ?? ErrorColor;
}

export const ConfalPyramidsColorThemeParams = {};
export type ConfalPyramidsColorThemeParams = typeof ConfalPyramidsColorThemeParams
export function getConfalPyramidsColorThemeParams(ctx: ThemeDataContext) {
    return ConfalPyramidsColorThemeParams; // TODO return copy
}

export function ConfalPyramidsColorTheme(ctx: ThemeDataContext, props: PD.Values<ConfalPyramidsColorThemeParams>): ColorTheme<ConfalPyramidsColorThemeParams> {
    function color(location: Location, isSecondary: boolean): Color {
        if (CPT.isLocation(location)) {
            const { pyramid, isLower } = location.data;
            return getConformerColor(pyramid.NtC, isLower);
        }

        return DefaultColor;
    }

    return {
        factory: ConfalPyramidsColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: TableLegend(iterableToArray(ColorMapping.entries()).map(([conformer, color]) => {
            return [conformer, color] as [string, Color];
        }).concat([
            [ 'Error', ErrorColor ],
            [ 'Unknown', DefaultColor ]
        ]))
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
        detach: (data) => data.structure && data.structure.models[0].customProperties.reference(ConfalPyramidsProvider.descriptor, false)
    }
};
