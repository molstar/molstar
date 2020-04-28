/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ShrakeRupleyComputationParams, AccessibleSurfaceArea } from './accessible-surface-area/shrake-rupley';
import { Structure, Unit } from '../../mol-model/structure';
import { CustomStructureProperty } from '../common/custom-structure-property';
import { CustomProperty } from '../common/custom-property';
import { QuerySymbolRuntime } from '../../mol-script/runtime/query/compiler';
import { CustomPropSymbol } from '../../mol-script/language/symbol';
import Type from '../../mol-script/language/type';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';

export const AccessibleSurfaceAreaParams = {
    ...ShrakeRupleyComputationParams
};
export type AccessibleSurfaceAreaParams = typeof AccessibleSurfaceAreaParams
export type AccessibleSurfaceAreaProps = PD.Values<AccessibleSurfaceAreaParams>

export const AccessibleSurfaceAreaSymbols = {
    isBuried: QuerySymbolRuntime.Dynamic(CustomPropSymbol('computed', 'accessible-surface-area.is-buried', Type.Bool),
        ctx => {
            if (!Unit.isAtomic(ctx.element.unit)) return false;
            const accessibleSurfaceArea = AccessibleSurfaceAreaProvider.get(ctx.element.structure).value;
            if (!accessibleSurfaceArea) return false;
            return AccessibleSurfaceArea.getFlag(ctx.element, accessibleSurfaceArea) === AccessibleSurfaceArea.Flag.Buried;
        }
    ),
    isAccessible: QuerySymbolRuntime.Dynamic(CustomPropSymbol('computed', 'accessible-surface-area.is-accessible', Type.Bool),
        ctx => {
            if (!Unit.isAtomic(ctx.element.unit)) return false;
            const accessibleSurfaceArea = AccessibleSurfaceAreaProvider.get(ctx.element.structure).value;
            if (!accessibleSurfaceArea) return false;
            return AccessibleSurfaceArea.getFlag(ctx.element, accessibleSurfaceArea) === AccessibleSurfaceArea.Flag.Accessible;
        }
    ),
};

export type AccessibleSurfaceAreaValue = AccessibleSurfaceArea

export const AccessibleSurfaceAreaProvider: CustomStructureProperty.Provider<AccessibleSurfaceAreaParams, AccessibleSurfaceAreaValue> = CustomStructureProperty.createProvider({
    label: 'Accessible Surface Area',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_accessible_surface_area',
        symbols: AccessibleSurfaceAreaSymbols,
        // TODO `cifExport`
    }),
    type: 'root',
    defaultParams: AccessibleSurfaceAreaParams,
    getParams: (data: Structure) => AccessibleSurfaceAreaParams,
    isApplicable: (data: Structure) => true,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<AccessibleSurfaceAreaProps>) => {
        const p = { ...PD.getDefaultValues(AccessibleSurfaceAreaParams), ...props };
        return { value: await AccessibleSurfaceArea.compute(data, p).runInContext(ctx.runtime) };
    }
});