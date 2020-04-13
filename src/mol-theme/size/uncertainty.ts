/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Unit, Bond, ElementIndex } from '../../mol-model/structure';
import { Location } from '../../mol-model/location';
import { SizeTheme } from '../size';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';

const Description = `Assigns a size reflecting the uncertainty or disorder of an element's position, e.g. B-factor or RMSF, depending on the data availability and experimental technique.`;

export const UncertaintySizeThemeParams = {
    bfactorFactor: PD.Numeric(0.1, { min: 0, max: 1, step: 0.01 }),
    rmsfFactor: PD.Numeric(0.05, { min: 0, max: 1, step: 0.01 }),
    baseSize: PD.Numeric(0.2, { min: 0, max: 2, step: 0.1 }),
};
export type UncertaintySizeThemeParams = typeof UncertaintySizeThemeParams
export function getUncertaintySizeThemeParams(ctx: ThemeDataContext) {
    return UncertaintySizeThemeParams; // TODO return copy
}

export function getUncertainty(unit: Unit, element: ElementIndex, props: PD.Values<UncertaintySizeThemeParams>): number {
    if (Unit.isAtomic(unit)) {
        return unit.model.atomicConformation.B_iso_or_equiv.value(element) * props.bfactorFactor;
    } else if (Unit.isSpheres(unit)) {
        return unit.model.coarseConformation.spheres.rmsf[element] * props.rmsfFactor;
    } else {
        return 0;
    }
}

export function UncertaintySizeTheme(ctx: ThemeDataContext, props: PD.Values<UncertaintySizeThemeParams>): SizeTheme<UncertaintySizeThemeParams> {
    function size(location: Location): number {
        let size = props.baseSize;
        if (StructureElement.Location.is(location)) {
            size += getUncertainty(location.unit, location.element, props);
        } else if (Bond.isLocation(location)) {
            size += getUncertainty(location.aUnit, location.aUnit.elements[location.aIndex], props);
        }
        return size;
    }

    return {
        factory: UncertaintySizeTheme,
        granularity: 'group',
        size,
        props,
        description: Description
    };
}

export const UncertaintySizeThemeProvider: SizeTheme.Provider<UncertaintySizeThemeParams, 'uncertainty'> = {
    name: 'uncertainty',
    label: 'Uncertainty/Disorder',
    category: '',
    factory: UncertaintySizeTheme,
    getParams: getUncertaintySizeThemeParams,
    defaultValues: PD.getDefaultValues(UncertaintySizeThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.models.some(m => m.atomicConformation.B_iso_or_equiv.isDefined || m.coarseHierarchy.isDefined)
};