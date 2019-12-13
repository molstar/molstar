/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomPropertyDescriptor, Structure } from '../../mol-model/structure';
import { RuntimeContext } from '../../mol-task';
import { DSSPComputationParams, DSSPComputationProps, computeUnitDSSP } from './secondary-structure/dssp';
import { SecondaryStructure } from '../../mol-model/structure/model/properties/seconday-structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Unit } from '../../mol-model/structure/structure';
import { CustomStructureProperty } from '../common/custom-property-registry';

// TODO get default params based on structure
// /**
//  * Attaches ComputedSecondaryStructure property when unavailable in sourceData and
//  * when not an archival file (i.e. no database_2.database_id field)
//  */
// export async function ensureSecondaryStructure(s: Structure) {
//     if (s.models.length === 1 && s.model && s.model.sourceData.kind === 'mmCIF') {
//         if (!s.model.sourceData.data.struct_conf.id.isDefined && !s.model.sourceData.data.struct_sheet_range.id.isDefined &&
//             !s.model.sourceData.data.database_2.database_id.isDefined
//         ) {
//             await ComputedSecondaryStructure.attach(s)
//         }
//     }
// }

export const SecondaryStructureParams = {
    type: PD.MappedStatic('mmcif', {
        'mmcif': PD.EmptyGroup({ label: 'mmCIF' }),
        'dssp': PD.Group(DSSPComputationParams, { label: 'DSSP', isFlat: true })
    }, { options: [['mmcif', 'mmCIF'], ['dssp', 'DSSP']] })
}
export type SecondaryStructureParams = typeof SecondaryStructureParams
export type SecondaryStructureProps = PD.Values<SecondaryStructureParams>

export type SecondaryStructureValue = Map<number, SecondaryStructure>

export const SecondaryStructureProvider: CustomStructureProperty.Provider<SecondaryStructureParams, SecondaryStructureValue> = CustomStructureProperty.createProvider({
    label: 'Secondary Structure',
    descriptor: CustomPropertyDescriptor({
        isStatic: true,
        name: 'molstar_computed_secondary_structure',
        // TODO `cifExport` and `symbol`
    }),
    defaultParams: SecondaryStructureParams,
    getParams: (data: Structure) => SecondaryStructureParams,
    isApplicable: (data: Structure) => true,
    compute: async (ctx: RuntimeContext, data: Structure, props: Partial<SecondaryStructureProps>) => {
        const p = { ...PD.getDefaultValues(SecondaryStructureParams), ...props }
        switch (p.type.name) {
            case 'dssp': return await computeDssp(data, p.type.params)
            case 'mmcif': return await computeMmcif(data)
        }
    }
})

async function computeDssp(structure: Structure, props: DSSPComputationProps): Promise<SecondaryStructureValue> {
    // TODO take inter-unit hbonds into account for bridge, ladder, sheet assignment
    // TODO store unit-only secStruc as custom unit property???
    // TODO use Zhang-Skolnik for CA alpha only parts or for coarse parts with per-residue elements
    const map = new Map<number, SecondaryStructure>()
    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const u = structure.unitSymmetryGroups[i].units[0]
        if (Unit.isAtomic(u)) {
            const secondaryStructure = await computeUnitDSSP(u, props)
            map.set(u.invariantId, secondaryStructure)
        }
    }
    return map
}

async function computeMmcif(structure: Structure): Promise<SecondaryStructureValue> {
    const map = new Map<number, SecondaryStructure>()
    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const u = structure.unitSymmetryGroups[i].units[0]
        if (Unit.isAtomic(u)) {
            map.set(u.invariantId, u.model.properties.secondaryStructure)
        }
    }
    return map
}