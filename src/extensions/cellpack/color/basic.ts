/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../../mol-util/color';
import { StructureElement, Unit, ElementIndex, Bond } from '../../../mol-model/structure';
import { Location } from '../../../mol-model/location';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { ColorNames } from '../../../mol-util/color/names';
import { EntityIdColorTheme, EntityIdColorThemeParams } from '../../../mol-theme/color/entity-id';

const DefaultCellPackColor = Color(0xEEEE00);
const Description = 'CellPack coloring.';

export const CellPackColorThemeParams = {
    ...EntityIdColorThemeParams,
};
export type CellPackColorThemeParams = typeof CellPackColorThemeParams
export function getCellPackColorThemeParams(ctx: ThemeDataContext) {
    return CellPackColorThemeParams; // TODO return copy
}

function getAtomicCompId(unit: Unit.Atomic, element: ElementIndex) {
    return unit.model.atomicHierarchy.atoms.label_comp_id.value(element);
}

const lipids = ['LIP'];

export function CellPackColorTheme(ctx: ThemeDataContext, props: PD.Values<CellPackColorThemeParams>): ColorTheme<CellPackColorThemeParams> {
    let color: LocationColor;

    if (ctx.structure) {
        const entityColor = EntityIdColorTheme(ctx, props).color;

        color = (location: Location, isSecondary: boolean): Color => {
            if (StructureElement.Location.is(location)) {
                if (Unit.isAtomic(location.unit)) {
                    const compId = getAtomicCompId(location.unit, location.element);
                    if (lipids.includes(compId)) {
                        return ColorNames.lightgrey;
                    }
                }
            } else if (Bond.isLocation(location)) {
                if (Unit.isAtomic(location.aUnit)) {
                    const compId = getAtomicCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
                    if (lipids.includes(compId)) {
                        return ColorNames.lightgrey;
                    }
                }
            }
            return entityColor(location, isSecondary);
        };
    } else {
        color = () => DefaultCellPackColor;
    }

    return {
        factory: CellPackColorTheme,
        granularity: 'group',
        preferSmoothing: false,
        color,
        props,
        description: Description,
    };
}

export const CellPackColorThemeProvider: ColorTheme.Provider<CellPackColorThemeParams, 'cellpack'> = {
    name: 'cellpack',
    label: 'CellPack',
    category: ColorTheme.Category.Misc,
    factory: CellPackColorTheme,
    getParams: getCellPackColorThemeParams,
    defaultValues: PD.getDefaultValues(CellPackColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
};