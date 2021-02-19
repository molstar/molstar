/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureProperties, StructureElement, Bond, Structure } from '../../mol-model/structure';
import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';

const DefaultList = 'many-distinct';
const DefaultColor = Color(0xFAFAFA);
const Description = 'Gives every chain a color based on its `label_entity_id` value.';

export const EntityIdColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type EntityIdColorThemeParams = typeof EntityIdColorThemeParams
export function getEntityIdColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(EntityIdColorThemeParams);
    return params;
}

function key(entityId: string, modelIndex: number) {
    return `${entityId}|${modelIndex}`;
}

function getEntityIdSerialMap(structure: Structure) {
    const map = new Map<string, number>();
    for (let i = 0, il = structure.models.length; i < il; ++i) {
        const { label_entity_id } = structure.models[i].atomicHierarchy.chains;
        for (let j = 0, jl = label_entity_id.rowCount; j < jl; ++j) {
            const k = key(label_entity_id.value(j), i);
            if (!map.has(k)) map.set(k, map.size);
        }
    }
    return map;
}

export function EntityIdColorTheme(ctx: ThemeDataContext, props: PD.Values<EntityIdColorThemeParams>): ColorTheme<EntityIdColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    if (ctx.structure) {
        const l = StructureElement.Location.create(ctx.structure.root);
        const entityIdSerialMap = getEntityIdSerialMap(ctx.structure.root);

        const labelTable = Array.from(entityIdSerialMap.keys());
        const valueLabel = (i: number) => labelTable[i];

        const palette = getPalette(entityIdSerialMap.size, props, { valueLabel });
        legend = palette.legend;

        color = (location: Location): Color => {
            let serial: number | undefined = undefined;
            if (StructureElement.Location.is(location)) {
                const atomId = StructureProperties.chain.label_entity_id(location);
                const modelIndex = location.structure.models.indexOf(location.unit.model);
                const k = key(atomId, modelIndex);
                serial = entityIdSerialMap.get(k);
            } else if (Bond.isLocation(location)) {
                l.unit = location.aUnit;
                l.element = location.aUnit.elements[location.aIndex];
                const atomId = StructureProperties.chain.label_entity_id(l);
                const modelIndex = l.structure.models.indexOf(l.unit.model);
                const k = key(atomId, modelIndex);
                serial = entityIdSerialMap.get(k);
            }
            return serial === undefined ? DefaultColor : palette.color(serial);
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: EntityIdColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend
    };
}

export const EntityIdColorThemeProvider: ColorTheme.Provider<EntityIdColorThemeParams, 'entity-id'> = {
    name: 'entity-id',
    label: 'Entity Id',
    category: ColorTheme.Category.Chain,
    factory: EntityIdColorTheme,
    getParams: getEntityIdColorThemeParams,
    defaultValues: PD.getDefaultValues(EntityIdColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};