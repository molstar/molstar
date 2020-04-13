/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureProperties, StructureElement, Bond, Structure } from '../../mol-model/structure';
import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';
import { Segmentation } from '../../mol-data/int';
import { ColorLists, getColorListFromName } from '../../mol-util/color/lists';

const DefaultList = 'dark-2';
const DefaultColor = Color(0xFAFAFA);
const Description = 'Gives every chain a color based on its `asym_id` value.';

export const ChainIdColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type ChainIdColorThemeParams = typeof ChainIdColorThemeParams
export function getChainIdColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(ChainIdColorThemeParams);
    if (ctx.structure) {
        if (getAsymIdSerialMap(ctx.structure.root).size > ColorLists[DefaultList].list.length) {
            params.palette.defaultValue.name = 'colors';
            params.palette.defaultValue.params = {
                ...params.palette.defaultValue.params,
                list: { kind: 'interpolate', colors: getColorListFromName(DefaultList).list }
            };
        }
    }
    return params;
}

function getAsymId(unit: Unit): StructureElement.Property<string> {
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            return StructureProperties.chain.label_asym_id;
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return StructureProperties.coarse.asym_id;
    }
}

function getAsymIdSerialMap(structure: Structure) {
    const map = new Map<string, number>();
    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const unit = structure.unitSymmetryGroups[i].units[0];
        const { model } = unit;
        if (Unit.isAtomic(unit)) {
            const { chainAtomSegments, chains } = model.atomicHierarchy;
            const chainIt = Segmentation.transientSegments(chainAtomSegments, unit.elements);
            while (chainIt.hasNext) {
                const { index: chainIndex } = chainIt.move();
                const asymId = chains.label_asym_id.value(chainIndex);
                if (!map.has(asymId)) map.set(asymId, map.size);
            }
        } else if (Unit.isCoarse(unit)) {
            const { chainElementSegments, asym_id } = Unit.isSpheres(unit)
                ? model.coarseHierarchy.spheres
                : model.coarseHierarchy.gaussians;
            const chainIt = Segmentation.transientSegments(chainElementSegments, unit.elements);
            while (chainIt.hasNext) {
                const { index: chainIndex } = chainIt.move();
                const elementIndex = chainElementSegments.offsets[chainIndex];
                const asymId = asym_id.value(elementIndex);
                if (!map.has(asymId)) map.set(asymId, map.size);
            }
        }
    }
    return map;
}

export function ChainIdColorTheme(ctx: ThemeDataContext, props: PD.Values<ChainIdColorThemeParams>): ColorTheme<ChainIdColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    if (ctx.structure) {
        const l = StructureElement.Location.create(ctx.structure);
        const asymIdSerialMap = getAsymIdSerialMap(ctx.structure.root);

        const labelTable = Array.from(asymIdSerialMap.keys());
        props.palette.params.valueLabel = (i: number) => labelTable[i];

        const palette = getPalette(asymIdSerialMap.size, props);
        legend = palette.legend;

        color = (location: Location): Color => {
            let serial: number | undefined = undefined;
            if (StructureElement.Location.is(location)) {
                const asym_id = getAsymId(location.unit);
                serial = asymIdSerialMap.get(asym_id(location));
            } else if (Bond.isLocation(location)) {
                const asym_id = getAsymId(location.aUnit);
                l.unit = location.aUnit;
                l.element = location.aUnit.elements[location.aIndex];
                serial = asymIdSerialMap.get(asym_id(l));
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