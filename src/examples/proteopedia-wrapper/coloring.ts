/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */


import { Unit, StructureProperties, StructureElement, Bond } from '../../mol-model/structure';

import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { ColorTheme, LocationColor } from '../../mol-theme/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { Column } from '../../mol-data/db';

const Description = 'Gives every chain a color from a list based on its `asym_id` value.';

export function createProteopediaCustomTheme(colors: number[]) {
    const ProteopediaCustomColorThemeParams = {
        colors: PD.ObjectList({ color: PD.Color(Color(0xffffff)) }, ({ color }) => Color.toHexString(color),
            { defaultValue: colors.map(c => ({ color: Color(c) })) })
    };
    type ProteopediaCustomColorThemeParams = typeof ProteopediaCustomColorThemeParams
    function getChainIdColorThemeParams(ctx: ThemeDataContext) {
        return ProteopediaCustomColorThemeParams; // TODO return copy
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

    function addAsymIds(map: Map<string, number>, data: Column<string>) {
        let j = map.size;
        for (let o = 0, ol = data.rowCount; o < ol; ++o) {
            const k = data.value(o);
            if (!map.has(k)) {
                map.set(k, j);
                j += 1;
            }
        }
    }

    function ProteopediaCustomColorTheme(ctx: ThemeDataContext, props: PD.Values<ProteopediaCustomColorThemeParams>): ColorTheme<ProteopediaCustomColorThemeParams> {
        let color: LocationColor;

        const colors = props.colors, colorCount = colors.length, defaultColor = colors[0].color;

        if (ctx.structure) {
            const l = StructureElement.Location.create(ctx.structure);
            const { models } = ctx.structure;
            const asymIdSerialMap = new Map<string, number>();
            for (let i = 0, il = models.length; i < il; ++i) {
                const m = models[i];
                addAsymIds(asymIdSerialMap, m.atomicHierarchy.chains.label_asym_id);
                if (m.coarseHierarchy.isDefined) {
                    addAsymIds(asymIdSerialMap, m.coarseHierarchy.spheres.asym_id);
                    addAsymIds(asymIdSerialMap, m.coarseHierarchy.gaussians.asym_id);
                }
            }

            color = (location: Location): Color => {
                if (StructureElement.Location.is(location)) {
                    const asym_id = getAsymId(location.unit);
                    const o = asymIdSerialMap.get(asym_id(location)) || 0;
                    return colors[o % colorCount].color;
                } else if (Bond.isLocation(location)) {
                    const asym_id = getAsymId(location.aUnit);
                    l.unit = location.aUnit;
                    l.element = location.aUnit.elements[location.aIndex];
                    const o = asymIdSerialMap.get(asym_id(l)) || 0;
                    return colors[o % colorCount].color;
                }
                return defaultColor;
            };
        } else {
            color = () => defaultColor;
        }

        return {
            factory: ProteopediaCustomColorTheme,
            granularity: 'group',
            color,
            props,
            description: Description,
            legend: undefined
        };
    }

    return {
        name: 'proteopedia-custom',
        label: 'Proteopedia Custom',
        category: 'Custom',
        factory: ProteopediaCustomColorTheme,
        getParams: getChainIdColorThemeParams,
        defaultValues: PD.getDefaultValues(ProteopediaCustomColorThemeParams),
        isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
    };
}