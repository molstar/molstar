/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssemblySymmetryQuery, AssemblySymmetryQueryVariables } from './graphql/types';
import query from './graphql/symmetry.gql';

import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { CustomPropertyDescriptor, Structure, Model } from '../../mol-model/structure';
import { Database as _Database, Column } from '../../mol-data/db'
import { GraphQLClient } from '../../mol-util/graphql-client';
import { CustomProperty } from '../common/custom-property';
import { NonNullableArray } from '../../mol-util/type-helpers';
import { CustomStructureProperty } from '../common/custom-structure-property';
import { MmcifFormat } from '../../mol-model-formats/structure/mmcif';
import { ReadonlyVec3 } from '../../mol-math/linear-algebra/3d/vec3';

const BiologicalAssemblyNames = new Set([
    'author_and_software_defined_assembly',
    'author_defined_assembly',
    'complete icosahedral assembly',
    'complete point assembly',
    'representative helical assembly',
    'software_defined_assembly'
])

export function isBiologicalAssembly(structure: Structure): boolean {
    if (!MmcifFormat.is(structure.models[0].sourceData)) return false
    const mmcif = structure.models[0].sourceData.data.db
    if (!mmcif.pdbx_struct_assembly.details.isDefined) return false
    const id = structure.units[0].conformation.operator.assembly?.id || ''
    if (id === '' || id === 'deposited') return true
    const indices = Column.indicesOf(mmcif.pdbx_struct_assembly.id, e => e === id)
    if (indices.length !== 1) return false
    const details = mmcif.pdbx_struct_assembly.details.value(indices[0])
    return BiologicalAssemblyNames.has(details)
}

export namespace AssemblySymmetry {
    export enum Tag {
        Cluster = 'rcsb-assembly-symmetry-cluster',
        Representation = 'rcsb-assembly-symmetry-3d'
    }

    export const DefaultServerUrl = 'https://data-beta.rcsb.org/graphql'

    export function isApplicable(structure?: Structure): boolean {
        return (
            !!structure && structure.models.length === 1 &&
            Model.isFromPdbArchive(structure.models[0]) &&
            isBiologicalAssembly(structure)
        )
    }

    export async function fetch(ctx: CustomProperty.Context, structure: Structure, props: AssemblySymmetryDataProps): Promise<AssemblySymmetryDataValue> {
        if (!isApplicable(structure)) return []

        const client = new GraphQLClient(props.serverUrl, ctx.fetch)
        const variables: AssemblySymmetryQueryVariables = {
            assembly_id: structure.units[0].conformation.operator.assembly?.id || 'deposited',
            entry_id: structure.units[0].model.entryId
        }
        const result = await client.request<AssemblySymmetryQuery>(ctx.runtime, query, variables)

        if (!result.assembly?.rcsb_struct_symmetry) {
            console.error('expected `rcsb_struct_symmetry` field')
            return []
        }
        return result.assembly.rcsb_struct_symmetry as AssemblySymmetryDataValue
    }

    /** Returns the index of the first non C1 symmetry or -1 */
    export function firstNonC1(assemblySymmetryData: AssemblySymmetryDataValue) {
        for (let i = 0, il = assemblySymmetryData.length; i < il; ++i) {
            if (assemblySymmetryData[i].symbol !== 'C1') return i
        }
        return -1
    }

    export type RotationAxes = ReadonlyArray<{ order: number, start: ReadonlyVec3, end: ReadonlyVec3 }>
    export function isRotationAxes(x: AssemblySymmetryValue['rotation_axes']): x is RotationAxes {
        return !!x && x.length > 0
    }
}

export function getSymmetrySelectParam(structure?: Structure) {
    const param = PD.Select<number>(0, [[0, 'First Symmetry']])
    if (structure) {
        const assemblySymmetryData = AssemblySymmetryDataProvider.get(structure).value
        if (assemblySymmetryData) {
            const options: [number, string][] = []
            for (let i = 0, il = assemblySymmetryData.length; i < il; ++i) {
                const { symbol, kind } = assemblySymmetryData[i]
                if (symbol !== 'C1') {
                    options.push([ i, `${i + 1}: ${symbol} ${kind}` ])
                }
            }
            if (options.length) {
                param.options = options
                param.defaultValue = options[0][0]
            }
        }
    }
    return param
}

//

export const AssemblySymmetryDataParams = {
    serverUrl: PD.Text(AssemblySymmetry.DefaultServerUrl, { description: 'GraphQL endpoint URL' })
}
export type AssemblySymmetryDataParams = typeof AssemblySymmetryDataParams
export type AssemblySymmetryDataProps = PD.Values<AssemblySymmetryDataParams>

export type AssemblySymmetryDataValue = NonNullableArray<NonNullable<NonNullable<AssemblySymmetryQuery['assembly']>['rcsb_struct_symmetry']>>

export const AssemblySymmetryDataProvider: CustomStructureProperty.Provider<AssemblySymmetryDataParams, AssemblySymmetryDataValue> = CustomStructureProperty.createProvider({
    label: 'Assembly Symmetry Data',
    descriptor: CustomPropertyDescriptor({
        name: 'rcsb_struct_symmetry_data',
        // TODO `cifExport` and `symbol`
    }),
    type: 'root',
    defaultParams: AssemblySymmetryDataParams,
    getParams: (data: Structure) => AssemblySymmetryDataParams,
    isApplicable: (data: Structure) => AssemblySymmetry.isApplicable(data),
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<AssemblySymmetryDataProps>) => {
        const p = { ...PD.getDefaultValues(AssemblySymmetryDataParams), ...props }
        return await AssemblySymmetry.fetch(ctx, data, p)
    }
})

//

function getAssemblySymmetryParams(data?: Structure) {
    return {
        ... AssemblySymmetryDataParams,
        symmetryIndex: getSymmetrySelectParam(data)
    }
}

export const AssemblySymmetryParams = getAssemblySymmetryParams()
export type AssemblySymmetryParams = typeof AssemblySymmetryParams
export type AssemblySymmetryProps = PD.Values<AssemblySymmetryParams>

export type AssemblySymmetryValue = AssemblySymmetryDataValue[0]

export const AssemblySymmetryProvider: CustomStructureProperty.Provider<AssemblySymmetryParams, AssemblySymmetryValue> = CustomStructureProperty.createProvider({
    label: 'Assembly Symmetry',
    descriptor: CustomPropertyDescriptor({
        name: 'rcsb_struct_symmetry',
        // TODO `cifExport` and `symbol`
    }),
    type: 'root',
    defaultParams: AssemblySymmetryParams,
    getParams: getAssemblySymmetryParams,
    isApplicable: (data: Structure) => AssemblySymmetry.isApplicable(data),
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<AssemblySymmetryProps>) => {
        const p = { ...PD.getDefaultValues(getAssemblySymmetryParams(data)), ...props }
        await AssemblySymmetryDataProvider.attach(ctx, data, p)
        const assemblySymmetryData = AssemblySymmetryDataProvider.get(data).value
        const assemblySymmetry = assemblySymmetryData?.[p.symmetryIndex]
        if (!assemblySymmetry) new Error(`No assembly symmetry found for index ${p.symmetryIndex}`)
        return assemblySymmetry
    }
})