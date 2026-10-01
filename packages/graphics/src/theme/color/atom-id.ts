/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureProperties, StructureElement, Bond, type Structure } from '@molstar/model/model/structure';
import { Color } from '@molstar/core/util/color';
import type { Location } from '@molstar/model/model/location';
import type { ColorTheme, LocationColor } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { getPaletteParams, getPalette } from '@molstar/core/util/color/palette';
import type { TableLegend, ScaleLegend } from '@molstar/core/util/legend';
import { ColorThemeCategory } from './categories.js';

const DefaultList = 'many-distinct';
const DefaultColor = Color(0xFAFAFA);
const Description = 'Gives every atom a color based on its `label_atom_id` value.';

export const AtomIdColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type AtomIdColorThemeParams = typeof AtomIdColorThemeParams
export function getAtomIdColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(AtomIdColorThemeParams);
    return params;
}

function getAtomIdSerialMap(structure: Structure) {
    const map = new Map<string, number>();
    for (const m of structure.models) {
        const { label_atom_id } = m.atomicHierarchy.atoms;
        for (let i = 0, il = label_atom_id.rowCount; i < il; ++i) {
            const id = label_atom_id.value(i);
            if (!map.has(id)) map.set(id, map.size);
        }
    }
    return map;
}

export function AtomIdColorTheme(ctx: ThemeDataContext, props: PD.Values<AtomIdColorThemeParams>): ColorTheme<AtomIdColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    if (ctx.structure) {
        const l = StructureElement.Location.create(ctx.structure.root);
        const atomIdSerialMap = getAtomIdSerialMap(ctx.structure.root);

        const labelTable = Array.from(atomIdSerialMap.keys());
        const valueLabel = (i: number) => labelTable[i];

        const palette = getPalette(atomIdSerialMap.size, props, { valueLabel });
        legend = palette.legend;

        color = (location: Location): Color => {
            let serial: number | undefined = undefined;
            if (StructureElement.Location.is(location)) {
                const id = StructureProperties.atom.label_atom_id(location);
                serial = atomIdSerialMap.get(id);
            } else if (Bond.isLocation(location)) {
                l.unit = location.aUnit;
                l.element = location.aUnit.elements[location.aIndex];
                const id = StructureProperties.atom.label_atom_id(l);
                serial = atomIdSerialMap.get(id);
            }
            return serial === undefined ? DefaultColor : palette.color(serial);
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: AtomIdColorTheme,
        granularity: 'group',
        preferSmoothing: true,
        color,
        props,
        description: Description,
        legend
    };
}

export const AtomIdColorThemeProvider: ColorTheme.Provider<AtomIdColorThemeParams, 'atom-id'> = {
    name: 'atom-id',
    label: 'Atom Id',
    category: ColorThemeCategory.Atom,
    factory: AtomIdColorTheme,
    getParams: getAtomIdColorThemeParams,
    defaultValues: PD.getDefaultValues(AtomIdColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};