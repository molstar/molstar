/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssemblySymmetryQuery, AssemblySymmetryQueryVariables } from '../graphql/types';
import { symmetry_gql } from '../graphql/symmetry.gql';

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Model, StructureSelection, QueryContext } from '../../../mol-model/structure';
import { Database as _Database, Column } from '../../../mol-data/db';
import { GraphQLClient } from '../../../mol-util/graphql-client';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { NonNullableArray } from '../../../mol-util/type-helpers';
import { CustomStructureProperty } from '../../../mol-model-props/common/custom-structure-property';
import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { ReadonlyVec3 } from '../../../mol-math/linear-algebra/3d/vec3';
import { SetUtils } from '../../../mol-util/set';
import { MolScriptBuilder as MS } from '../../../mol-script/language/builder';
import { compile } from '../../../mol-script/runtime/query/compiler';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';
import { Asset } from '../../../mol-util/assets';

const BiologicalAssemblyNames = new Set([
    'author_and_software_defined_assembly',
    'author_defined_assembly',
    'complete icosahedral assembly',
    'complete point assembly',
    'representative helical assembly',
    'software_defined_assembly'
]);

export function isBiologicalAssembly(structure: Structure): boolean {
    if (!MmcifFormat.is(structure.models[0].sourceData)) return false;
    const mmcif = structure.models[0].sourceData.data.db;
    if (!mmcif.pdbx_struct_assembly.details.isDefined) return false;
    const id = structure.units[0].conformation.operator.assembly?.id || '';
    if (id === '') return true;
    const indices = Column.indicesOf(mmcif.pdbx_struct_assembly.id, e => e === id);
    if (indices.length !== 1) return false;
    const details = mmcif.pdbx_struct_assembly.details.value(indices[0]);
    return BiologicalAssemblyNames.has(details);
}

export namespace AssemblySymmetry {
    export enum Tag {
        Cluster = 'rcsb-assembly-symmetry-cluster',
        Representation = 'rcsb-assembly-symmetry-3d'
    }

    export const DefaultServerUrl = 'https://data.rcsb.org/graphql'; // Alternative: 'https://www.ebi.ac.uk/pdbe/aggregated-api/pdb/symmetry' (if serverType is 'pdbe')

    export function isApplicable(structure?: Structure): boolean {
        return (
            !!structure && structure.models.length === 1 &&
            Model.hasPdbId(structure.models[0]) &&
            isBiologicalAssembly(structure)
        );
    }

    export async function fetch(ctx: CustomProperty.Context, structure: Structure, props: AssemblySymmetryDataProps): Promise<CustomProperty.Data<AssemblySymmetryDataValue>> {
        if (!isApplicable(structure)) return { value: [] };

        if (props.serverType === 'pdbe') return fetch_PDBe(ctx, structure, props);

        const client = new GraphQLClient(props.serverUrl, ctx.assetManager);
        const variables: AssemblySymmetryQueryVariables = {
            assembly_id: structure.units[0].conformation.operator.assembly?.id || '',
            entry_id: structure.units[0].model.entryId
        };
        const result = await client.request(ctx.runtime, symmetry_gql, variables);
        let value: AssemblySymmetryDataValue = [];

        if (!result.data.assembly?.rcsb_struct_symmetry) {
            console.error('expected `rcsb_struct_symmetry` field');
        } else {
            value = result.data.assembly.rcsb_struct_symmetry as AssemblySymmetryDataValue;
        }
        return { value, assets: [result] };
    }

    async function fetch_PDBe(ctx: CustomProperty.Context, structure: Structure, props: AssemblySymmetryDataProps): Promise<CustomProperty.Data<AssemblySymmetryDataValue>> {
        const assembly_id = structure.units[0].conformation.operator.assembly?.id || '-1'; // should use '' instead of '-1' but the API does not support non-number assembly_id
        const entry_id = structure.units[0].model.entryId.toLowerCase();
        const url = `${props.serverUrl}/${entry_id}?assembly_id=${assembly_id}`;
        const asset = Asset.getUrlAsset(ctx.assetManager, url);
        let dataWrapper: Asset.Wrapper<'json'>;
        try {
            dataWrapper = await ctx.assetManager.resolve(asset, 'json').runInContext(ctx.runtime);
        } catch (err) {
            // PDBe API returns 404 when there are no symmetries -> treat as success with empty json in body
            if (`${err}`.includes('404')) { // dirrrty
                dataWrapper = Asset.Wrapper({}, asset, ctx.assetManager);
            } else {
                throw err;
            }
        }
        const data = dataWrapper.data;

        const value: AssemblySymmetryDataValue = (data[entry_id] ?? []).map((v: any) => ({
            kind: 'Global Symmetry',
            oligomeric_state: v.oligomeric_state,
            stoichiometry: [v.stoichiometry],
            symbol: v.symbol,
            type: v.type,
            clusters: [],
            rotation_axes: v.rotation_axes,
        }));

        return { value, assets: [dataWrapper] };
    }

    /** Returns the index of the first non C1 symmetry or -1 */
    export function firstNonC1(assemblySymmetryData: AssemblySymmetryDataValue) {
        for (let i = 0, il = assemblySymmetryData.length; i < il; ++i) {
            if (assemblySymmetryData[i].symbol !== 'C1') return i;
        }
        return -1;
    }

    export type RotationAxes = ReadonlyArray<{ order: number, start: ReadonlyVec3, end: ReadonlyVec3 }>
    export function isRotationAxes(x: AssemblySymmetryValue['rotation_axes']): x is RotationAxes {
        return !!x && x.length > 0;
    }

    export function getAsymIds(assemblySymmetry: AssemblySymmetryValue) {
        const asymIds = new Set<string>();
        for (const c of assemblySymmetry.clusters) {
            if (!c?.members) continue;
            for (const m of c.members) {
                if (m?.asym_id) asymIds.add(m.asym_id);
            }
        }
        return SetUtils.toArray(asymIds);
    }

    function getAsymIdsStructure(structure: Structure, asymIds: string[]) {
        const query = MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
                'chain-test': MS.core.set.has([MS.set(...asymIds), MS.ammp('label_asym_id')])
            })
        ]);
        const compiled = compile<StructureSelection>(query);
        const result = compiled(new QueryContext(structure));
        return StructureSelection.unionStructure(result);
    }

    /** Returns structure limited to all cluster member chains */
    export function getStructure(structure: Structure, assemblySymmetry: AssemblySymmetryValue) {
        const asymIds = AssemblySymmetry.getAsymIds(assemblySymmetry);
        return asymIds.length > 0 ? getAsymIdsStructure(structure, asymIds) : structure;
    }
}

export function getSymmetrySelectParam(structure?: Structure) {
    const param = PD.Select<number>(0, [[0, 'First Symmetry']]);
    if (structure) {
        const assemblySymmetryData = AssemblySymmetryDataProvider.get(structure).value;
        if (assemblySymmetryData) {
            const options: [number, string][] = [
                [-1, 'Off']
            ];
            for (let i = 0, il = assemblySymmetryData.length; i < il; ++i) {
                const { symbol, kind } = assemblySymmetryData[i];
                if (symbol !== 'C1') {
                    options.push([i, `${i + 1}: ${symbol} ${kind}`]);
                }
            }
            if (options.length > 1) {
                param.options = options;
                param.defaultValue = options[1][0];
            } else {
                options.length = 0;
            }
        }
    }
    return param;
}

//

export const AssemblySymmetryDataParams = {
    serverType: PD.Select('rcsb', [['rcsb', 'RCSB'], ['pdbe', 'PDBe']] as const),
    serverUrl: PD.Text(AssemblySymmetry.DefaultServerUrl, { description: 'GraphQL endpoint URL (if server type is RCSB) or PDBe API endpoint URL (if server type is PDBe)' })
};
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
        const p = { ...PD.getDefaultValues(AssemblySymmetryDataParams), ...props };
        return await AssemblySymmetry.fetch(ctx, data, p);
    }
});

//

function getAssemblySymmetryParams(data?: Structure) {
    return {
        ...AssemblySymmetryDataParams,
        symmetryIndex: getSymmetrySelectParam(data)
    };
}

export const AssemblySymmetryParams = getAssemblySymmetryParams();
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
        const p = { ...PD.getDefaultValues(getAssemblySymmetryParams(data)), ...props };
        await AssemblySymmetryDataProvider.attach(ctx, data, p);
        const assemblySymmetryData = AssemblySymmetryDataProvider.get(data).value;
        const assemblySymmetry = assemblySymmetryData?.[p.symmetryIndex];
        if (!assemblySymmetry) new Error(`No assembly symmetry found for index ${p.symmetryIndex}`);
        return { value: assemblySymmetry };
    }
});