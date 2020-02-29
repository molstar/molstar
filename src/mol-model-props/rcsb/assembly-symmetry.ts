/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssemblySymmetryQuery, AssemblySymmetryQueryVariables } from './graphql/types';
import query from './graphql/symmetry.gql';

import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { CustomPropertyDescriptor, Structure } from '../../mol-model/structure';
import { Database as _Database, Column } from '../../mol-data/db'
import { GraphQLClient } from '../../mol-util/graphql-client';
import { CustomProperty } from '../common/custom-property';
import { NonNullableArray } from '../../mol-util/type-helpers';
import { CustomStructureProperty } from '../common/custom-structure-property';
import { MmcifFormat } from '../../mol-model-formats/structure/mmcif';

const BiologicalAssemblyNames = new Set([
    'author_and_software_defined_assembly',
    'author_defined_assembly',
    'complete icosahedral assembly',
    'complete point assembly',
    'representative helical assembly',
    'software_defined_assembly'
])

export namespace AssemblySymmetry {
    export enum Tag {
        Cluster = 'rcsb-assembly-symmetry-cluster',
        Representation = 'rcsb-assembly-symmetry-3d'
    }

    export const DefaultServerUrl = 'https://data-beta.rcsb.org/graphql'

    export function isApplicable(structure?: Structure): boolean {
        // check if structure is from pdb entry
        if (!structure || structure.models.length !== 1 || !MmcifFormat.is(structure.models[0].sourceData) || (!structure.models[0].sourceData.data.db.database_2.database_id.isDefined &&
        structure.models[0].entryId.length !== 4)) return false

        // check if assembly is 'biological'
        const mmcif = structure.models[0].sourceData.data.db
        if (!mmcif.pdbx_struct_assembly.details.isDefined) return false
        const id = structure.units[0].conformation.operator.assembly.id
        if (id === '' || id === 'deposited') return true
        const indices = Column.indicesOf(mmcif.pdbx_struct_assembly.id, e => e === id)
        if (indices.length !== 1) return false
        const details = mmcif.pdbx_struct_assembly.details.value(indices[0])
        return BiologicalAssemblyNames.has(details)
    }

    export async function fetch(ctx: CustomProperty.Context, structure: Structure, props: AssemblySymmetryProps): Promise<AssemblySymmetryValue> {
        if (!isApplicable(structure)) return []

        const client = new GraphQLClient(props.serverUrl, ctx.fetch)
        const variables: AssemblySymmetryQueryVariables = {
            assembly_id: structure.units[0].conformation.operator.assembly.id || 'deposited',
            entry_id: structure.units[0].model.entryId
        }
        const result = await client.request<AssemblySymmetryQuery>(ctx.runtime, query, variables)

        if (!result.assembly?.rcsb_struct_symmetry) {
            console.error('expected `rcsb_struct_symmetry` field')
            return []
        }
        const symmetry = result.assembly.rcsb_struct_symmetry as AssemblySymmetryValue
        return symmetry.filter(s => s.symbol !== 'C1')
    }
}

export function getSymmetrySelectParam(structure?: Structure) {
    const param = PD.Select<number>(0, [[0, 'No Symmetries']])
    if (structure) {
        const assemblySymmetry = AssemblySymmetryProvider.get(structure).value
        if (assemblySymmetry) {
            const options: [number, string][] = []
            for (let i = 0, il = assemblySymmetry.length; i < il; ++i) {
                const { symbol, kind } = assemblySymmetry[i]
                options.push([ i, `${i + 1}: ${symbol} ${kind}` ])
            }
            if (options.length) {
                param.options = options
                param.defaultValue = options[0][0]
            }
        }
    }
    return param
}

export const AssemblySymmetryParams = {
    serverUrl: PD.Text(AssemblySymmetry.DefaultServerUrl, { description: 'GraphQL endpoint URL' })
}
export type AssemblySymmetryParams = typeof AssemblySymmetryParams
export type AssemblySymmetryProps = PD.Values<AssemblySymmetryParams>

export type AssemblySymmetryValue = NonNullableArray<NonNullable<NonNullable<AssemblySymmetryQuery['assembly']>['rcsb_struct_symmetry']>>

export const AssemblySymmetryProvider: CustomStructureProperty.Provider<AssemblySymmetryParams, AssemblySymmetryValue> = CustomStructureProperty.createProvider({
    label: 'Assembly Symmetry',
    descriptor: CustomPropertyDescriptor({
        name: 'rcsb_struct_symmetry',
        // TODO `cifExport` and `symbol`
    }),
    type: 'root',
    defaultParams: AssemblySymmetryParams,
    getParams: (data: Structure) => AssemblySymmetryParams,
    isApplicable: (data: Structure) => AssemblySymmetry.isApplicable(data),
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<AssemblySymmetryProps>) => {
        const p = { ...PD.getDefaultValues(AssemblySymmetryParams), ...props }
        return await AssemblySymmetry.fetch(ctx, data, p)
    }
})