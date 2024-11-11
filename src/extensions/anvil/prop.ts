/**
 * Copyright (c) 2020-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Structure, StructureProperties, Unit } from '../../mol-model/structure';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { ANVILParams, ANVILProps, computeANVIL, isInMembranePlane } from './algorithm';
import { CustomStructureProperty } from '../../mol-model-props/common/custom-structure-property';
import { CustomProperty } from '../../mol-model-props/common/custom-property';
import { Vec3 } from '../../mol-math/linear-algebra';
import { QuerySymbolRuntime } from '../../mol-script/runtime/query/base';
import { CustomPropSymbol } from '../../mol-script/language/symbol';
import { Type } from '../../mol-script/language/type';

export const MembraneOrientationParams = {
    ...ANVILParams
};
export type MembraneOrientationParams = typeof MembraneOrientationParams
export type MembraneOrientationProps = PD.Values<MembraneOrientationParams>

export { MembraneOrientation };

interface MembraneOrientation {
    // point in membrane boundary
    readonly planePoint1: Vec3,
    // point in opposite side of membrane boundary
    readonly planePoint2: Vec3,
    // normal vector of membrane layer
    readonly normalVector: Vec3,
    // the radius of the membrane layer
    readonly radius: number,
    readonly centroid: Vec3
}

namespace MembraneOrientation {
    export enum Tag {
        Representation = 'membrane-orientation-3d'
    }

    const pos = Vec3();
    export const symbols = {
        isTransmembrane: QuerySymbolRuntime.Dynamic(CustomPropSymbol('computed', 'membrane-orientation.is-transmembrane', Type.Bool),
            ctx => {
                const { unit, structure } = ctx.element;
                const { x, y, z } = StructureProperties.atom;
                if (!Unit.isAtomic(unit)) return 0;
                const membraneOrientation = MembraneOrientationProvider.get(structure).value;
                if (!membraneOrientation) return 0;
                Vec3.set(pos, x(ctx.element), y(ctx.element), z(ctx.element));
                const { normalVector, planePoint1, planePoint2 } = membraneOrientation!;
                return isInMembranePlane(pos, normalVector, planePoint1, planePoint2);
            })
    };
}

export const MembraneOrientationProvider: CustomStructureProperty.Provider<MembraneOrientationParams, MembraneOrientation> = CustomStructureProperty.createProvider({
    label: 'Membrane Orientation',
    descriptor: CustomPropertyDescriptor({
        name: 'anvil_computed_membrane_orientation',
        symbols: MembraneOrientation.symbols,
        // TODO `cifExport`
    }),
    type: 'root',
    defaultParams: MembraneOrientationParams,
    getParams: (data: Structure) => MembraneOrientationParams,
    isApplicable,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<MembraneOrientationProps>) => {
        const p = { ...PD.getDefaultValues(MembraneOrientationParams), ...props };
        try {
            return { value: await computeAnvil(ctx, data, p) };
        } catch (e) {
            // the "Residues Embedded in Membrane" symbol may bypass isApplicable() checks
            console.warn('Failed to predict membrane orientation. This happens for short peptides and entries without amino acids.');
            return { value: undefined };
        }
    }
});

function isApplicable(structure: Structure) {
    if (!structure.isAtomic) return false;

    for (const model of structure.models) {
        const { byEntityKey } = model.sequence;
        for (const key of Object.keys(byEntityKey)) {
            const { kind, length } = byEntityKey[+key].sequence;
            if (kind !== 'protein') continue; // can only process protein chains
            if (length >= 15) return true; // short peptides might fail
        }
    }
    return false;
}

async function computeAnvil(ctx: CustomProperty.Context, data: Structure, props: Partial<ANVILProps>): Promise<MembraneOrientation> {
    const p = { ...PD.getDefaultValues(ANVILParams), ...props };
    return await computeANVIL(data, p).runInContext(ctx.runtime);
}