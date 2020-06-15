/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureProperties, StructureElement, Bond, Structure, Model } from '../../mol-model/structure';
import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';

const DefaultList = 'many-distinct';
const DefaultColor = Color(0xFAFAFA);
const Description = 'Gives every chain a color based on its `asym_id` value.';

export const ChainIdColorThemeParams = {
    asymId: PD.Select('auth', PD.arrayToOptions<AsymIdType>(['auth', 'label'])),
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type ChainIdColorThemeParams = typeof ChainIdColorThemeParams
export function getChainIdColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(ChainIdColorThemeParams);
    return params;
}

type AsymIdType = 'auth' | 'label';

function getAsymId(unit: Unit, type: AsymIdType): StructureElement.Property<string> {
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            return type === 'auth'
                ? StructureProperties.chain.auth_asym_id
                : StructureProperties.chain.label_asym_id;
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return StructureProperties.coarse.asym_id;
    }
}

function getAsymIdKey(location: StructureElement.Location, type: AsymIdType) {
    const asymId = getAsymId(location.unit, type)(location);
    return location.structure.models.length > 1
        ? getKey(location.unit.model, asymId)
        : asymId;
}

function getKey(model: Model, asymId: string) {
    return `${asymId} | ${(Model.Index.get(model).value || 0) + 1}`;
}

function getAsymIdSerialMap(structure: Structure, type: AsymIdType) {
    const map = new Map<string, number>();
    for (const m of structure.models) {
        const asymIdOffset = Model.AsymIdOffset.get(m).value;
        const offset = (type === 'auth' ? asymIdOffset?.auth : asymIdOffset?.label) || 0;
        m.properties.structAsymMap.forEach(({ auth_id }, label_id) => {
            const asymId = type === 'auth' ? auth_id : label_id;
            const k = structure.models.length > 1
                ? getKey(m, asymId)
                : asymId;
            if (!map.has(k)) map.set(k, map.size + offset);
        });
    }
    return map;
}

export function ChainIdColorTheme(ctx: ThemeDataContext, props: PD.Values<ChainIdColorThemeParams>): ColorTheme<ChainIdColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    if (ctx.structure) {
        const l = StructureElement.Location.create(ctx.structure);
        const asymIdSerialMap = getAsymIdSerialMap(ctx.structure.root, props.asymId);

        const labelTable = Array.from(asymIdSerialMap.keys());
        props.palette.params.valueLabel = (i: number) => labelTable[i];

        const palette = getPalette(asymIdSerialMap.size, props);
        legend = palette.legend;

        color = (location: Location): Color => {
            let serial: number | undefined = undefined;
            if (StructureElement.Location.is(location)) {
                const k = getAsymIdKey(location, props.asymId);
                serial = asymIdSerialMap.get(k);
            } else if (Bond.isLocation(location)) {
                l.unit = location.aUnit;
                l.element = location.aUnit.elements[location.aIndex];
                const k = getAsymIdKey(l, props.asymId);
                serial = asymIdSerialMap.get(k);
            }
            return serial === undefined ? DefaultColor : palette.color(serial);
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: ChainIdColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend
    };
}

export const ChainIdColorThemeProvider: ColorTheme.Provider<ChainIdColorThemeParams, 'chain-id'> = {
    name: 'chain-id',
    label: 'Chain Id',
    category: ColorTheme.Category.Chain,
    factory: ChainIdColorTheme,
    getParams: getChainIdColorThemeParams,
    defaultValues: PD.getDefaultValues(ChainIdColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};