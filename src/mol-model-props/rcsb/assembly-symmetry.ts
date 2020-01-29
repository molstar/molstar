/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssemblySymmetryQuery, AssemblySymmetryQueryVariables } from './graphql/types';
import query from './graphql/symmetry.gql';

import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { CustomPropertyDescriptor, Structure } from '../../mol-model/structure';
import { Database as _Database } from '../../mol-data/db'
import { GraphQLClient } from '../../mol-util/graphql-client';
import { CustomProperty } from '../common/custom-property';
import { NonNullableArray } from '../../mol-util/type-helpers';
import { CustomStructureProperty } from '../common/custom-structure-property';

export namespace AssemblySymmetry {
    export const DefaultServerUrl = 'http://data-staging.rcsb.org/graphql'

    export function isApplicable(structure?: Structure): boolean {
        return (
            !!structure &&
            structure.models.length === 1 &&
            structure.models[0].sourceData.kind === 'mmCIF' &&
            (structure.models[0].sourceData.data.database_2.database_id.isDefined ||
                structure.models[0].entryId.length === 4)
        )
    }

    export async function fetch(ctx: CustomProperty.Context, structure: Structure, props: AssemblySymmetryProps): Promise<AssemblySymmetryValue> {
        const client = new GraphQLClient(props.serverUrl, ctx.fetch)
        const variables: AssemblySymmetryQueryVariables = {
            assembly_id: structure.units[0].conformation.operator.assembly.id,
            entry_id: structure.units[0].model.entryId
        }
        const result = await client.request<AssemblySymmetryQuery>(ctx.runtime, query, variables)

        if (!result.assembly?.rcsb_struct_symmetry) {
            throw new Error('missing fields')
        }
        return result.assembly.rcsb_struct_symmetry as AssemblySymmetryValue
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
        isStatic: true,
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