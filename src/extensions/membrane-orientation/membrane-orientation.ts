/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Structure, StructureProperties } from '../../mol-model/structure';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { ANVILParams, ANVILProps, computeANVIL, isInMembranePlane } from './ANVIL';
import { CustomStructureProperty } from '../../mol-model-props/common/custom-structure-property';
import { CustomProperty } from '../../mol-model-props/common/custom-property';
import { AccessibleSurfaceAreaProvider } from '../../mol-model-props/computed/accessible-surface-area';
import { Vec3 } from '../../mol-math/linear-algebra';
import { QuerySymbolRuntime } from '../../mol-script/runtime/query/base';
import { CustomPropSymbol } from '../../mol-script/language/symbol';
import Type from '../../mol-script/language/type';
import { StructureSelectionQuery, StructureSelectionCategory } from '../../mol-plugin-state/helpers/structure-selection-query';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';

export const MembraneOrientationParams = {
    ...ANVILParams
};
export type MembraneOrientationParams = typeof MembraneOrientationParams
export type MembraneOrientationProps = PD.Values<MembraneOrientationParams>

interface MembraneOrientation {
    // point in membrane boundary
    readonly p1: Vec3,
    // point in opposite side of membrane boundary
    readonly p2: Vec3,
    // normal vector of membrane layer
    readonly normal: Vec3,
    // the radius of the membrane layer
    readonly radius: number,
    readonly centroid: Vec3
}

namespace MembraneOrientation {
    const pos = Vec3();
    export const symbols = {
        isTransmembrane: QuerySymbolRuntime.Dynamic(CustomPropSymbol('membrane-orientation', 'is-transmembrane', Type.Bool),
            ctx => {
                const structure = ctx.currentStructure;
                const { x, y, z } = StructureProperties.atom;
                if (!structure.isAtomic) return false;
                const membraneOrientation = MembraneOrientationProvider.get(structure).value;
                if (!membraneOrientation) return false;
                Vec3.set(pos, x(ctx.element), y(ctx.element), z(ctx.element));
                const { normal, p1, p2 } = membraneOrientation!;
                return isInMembranePlane(pos, normal, p1, p2);
            })
    };

    export const isTransmembrane = StructureSelectionQuery('Residues embedded in membrane', MS.struct.modifier.union([
        MS.struct.modifier.wholeResidues([
            MS.struct.modifier.union([
                MS.struct.generator.atomGroups({
                    'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
                    'atom-test': symbols.isTransmembrane.symbol(),
                })
            ])
        ])
    ]), {
        description: 'Select residues that are embedded between the membrane layers.',
        category: StructureSelectionCategory.Residue,
        ensureCustomProperties: (ctx, structure) => {
            return MembraneOrientationProvider.attach(ctx, structure);
        }
    });
}

export { MembraneOrientation };

export const MembraneOrientationProvider: CustomStructureProperty.Provider<MembraneOrientationParams, MembraneOrientation> = CustomStructureProperty.createProvider({
    label: 'Membrane Orientation',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_computed_membrane_orientation',
        symbols: MembraneOrientation.symbols
        // TODO `cifExport`
    }),
    type: 'root',
    defaultParams: MembraneOrientationParams,
    getParams: (data: Structure) => MembraneOrientationParams,
    isApplicable: (data: Structure) => true,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<MembraneOrientationProps>) => {
        const p = { ...PD.getDefaultValues(MembraneOrientationParams), ...props };
        return { value: await computeAnvil(ctx, data, p) };
    }
});

async function computeAnvil(ctx: CustomProperty.Context, data: Structure, props: Partial<ANVILProps>): Promise<MembraneOrientation> {
    await AccessibleSurfaceAreaProvider.attach(ctx, data);
    const p = { ...PD.getDefaultValues(ANVILParams), ...props };
    return await computeANVIL(data, p).runInContext(ctx.runtime);
}