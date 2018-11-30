/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssemblySymmetry as AssemblySymmetryGraphQL } from './graphql/types';
import query from './graphql/symmetry.gql';

import { Model, ModelPropertyDescriptor } from 'mol-model/structure';
import { CifWriter } from 'mol-io/writer/cif';
import { Database as _Database, Column, Table } from 'mol-data/db'
import { Category } from 'mol-io/writer/cif/encoder';
import { Tensor } from 'mol-math/linear-algebra';
import { CifExportContext } from 'mol-model/structure/export/mmcif';
import { toTable } from 'mol-io/reader/cif/schema';
import { CifCategory } from 'mol-io/reader/cif';
import { PropertyWrapper } from 'mol-model-props/common/wrapper';
import { Task, RuntimeContext } from 'mol-task';
import { GraphQLClient } from 'mol-util/graphql-client';
import { ajaxGet } from 'mol-util/data-source';

const { str, int, float, Aliased, Vector, List } = Column.Schema;

function getInstance(name: keyof AssemblySymmetry.Schema): (ctx: CifExportContext) => CifWriter.Category.Instance<any, any> {
    return function(ctx: CifExportContext) {
        const assemblySymmetry = AssemblySymmetry.get(ctx.structures[0].model);
        return assemblySymmetry ? Category.ofTable(assemblySymmetry.db[name]) : CifWriter.Category.Empty;
    }
}

function getCategory(name: keyof AssemblySymmetry.Schema) {
    return { name, instance: getInstance(name) }
}

function createDatabaseFromJson(assemblies: ReadonlyArray<AssemblySymmetryGraphQL.Assemblies>): AssemblySymmetry.Database {
    const Schema = AssemblySymmetry.Schema

    const featureRows: Table.Row<typeof Schema.rcsb_assembly_symmetry>[] = []
    const clusterRows: Table.Row<typeof Schema.rcsb_assembly_symmetry_cluster>[] = []
    const clusterMemberRows: Table.Row<typeof Schema.rcsb_assembly_symmetry_cluster_member>[] = []
    const axisRows: Table.Row<typeof Schema.rcsb_assembly_symmetry_axis>[] = []

    let id = 1 // start feature ids at 1
    let clusterCount = 0
    for (let i = 0, il = assemblies.length; i < il; ++i) {
        const { pdbx_struct_assembly, rcsb_struct_symmetry, rcsb_struct_symmetry_provenance } = assemblies[i]
        if (!pdbx_struct_assembly || !rcsb_struct_symmetry ||!rcsb_struct_symmetry_provenance) continue
        const assembly_id = pdbx_struct_assembly.id
        for (let j = 0, jl = rcsb_struct_symmetry.length; j < jl; ++j) {
            const rss = rcsb_struct_symmetry[j]! // TODO upstream, array members should not be nullable
            featureRows.push({
                id,
                assembly_id,
                provenance: rcsb_struct_symmetry_provenance,
                type: rss.type,
                stoichiometry: rss.stoichiometry as string[],  // TODO upstream, array members should not be nullable
                kind: rss.kind,
                symbol: rss.symbol,
                oligomeric_state: rss.oligomeric_state
            })

            const clusters = rss.clusters
            if (clusters) {
                for (let k = 0, kl = clusters.length; k < kl; ++k) {
                    const c = clusters[k]! // TODO upstream, array members should not be nullable
                    const cluster_id = clusterCount + k
                    clusterRows.push({
                        id: cluster_id,
                        symmetry_id: id,
                        avg_rmsd: c.avg_rmsd || 0, // TODO upstream, should not be nullable, or???
                    })
                    for (let l = 0, ll = c.members.length; l < ll; ++l) {
                        const m = c.members[l]! // TODO upstream, array members should not be nullable
                        clusterMemberRows.push({
                            cluster_id: cluster_id,
                            asym_id: m.asym_id,
                            pdbx_struct_oper_list_ids: (m.pdbx_struct_oper_list_ids || []) as string[]
                        })
                    }
                }
                clusterCount += clusters.length
            }

            const axes = rss.rotation_axes
            if (axes) {
                for (let k = 0, kl = axes.length; k < kl; ++k) {
                    const a = axes[k]!
                    axisRows.push({
                        symmetry_id: id,
                        order: a.order!,  // TODO upstream, should not be nullable, or???
                        start: a.start as Tensor.Data, // TODO upstream, array members should not be nullable
                        end: a.end as Tensor.Data // TODO upstream, array members should not be nullable
                    })
                }
            }

            id += 1
        }
    }

    return _Database.ofTables('assembly_symmetry', Schema, {
        rcsb_assembly_symmetry: Table.ofRows(Schema.rcsb_assembly_symmetry, featureRows),
        rcsb_assembly_symmetry_cluster: Table.ofRows(Schema.rcsb_assembly_symmetry_cluster, clusterRows),
        rcsb_assembly_symmetry_cluster_member: Table.ofRows(Schema.rcsb_assembly_symmetry_cluster_member, clusterMemberRows),
        rcsb_assembly_symmetry_axis: Table.ofRows(Schema.rcsb_assembly_symmetry_axis, axisRows)
    })
}

function createDatabaseFromCif(model: Model): AssemblySymmetry.Database {
    const Schema = AssemblySymmetry.Schema

    const rcsb_assembly_symmetry = toTable(Schema.rcsb_assembly_symmetry, model.sourceData.frame.categories.rcsb_assembly_symmetry_feature)

    let rcsb_assembly_symmetry_cluster
    if (model.sourceData.frame.categoryNames.includes('rcsb_assembly_symmetry_cluster')) {
        rcsb_assembly_symmetry_cluster = toTable(Schema.rcsb_assembly_symmetry_cluster, model.sourceData.frame.categories.rcsb_assembly_symmetry_cluster)
    } else {
        rcsb_assembly_symmetry_cluster = toTable(Schema.rcsb_assembly_symmetry_cluster, CifCategory.empty as any)
    }

    let rcsb_assembly_symmetry_cluster_member
    if (model.sourceData.frame.categoryNames.includes('rcsb_assembly_symmetry_cluster_member')) {
        rcsb_assembly_symmetry_cluster_member = toTable(Schema.rcsb_assembly_symmetry_cluster_member, model.sourceData.frame.categories.rcsb_assembly_symmetry_cluster_member)
    } else {
        rcsb_assembly_symmetry_cluster_member = toTable(Schema.rcsb_assembly_symmetry_cluster_member, CifCategory.empty as any)
    }

    let rcsb_assembly_symmetry_axis
    if (model.sourceData.frame.categoryNames.includes('rcsb_assembly_symmetry_axis')) {
        rcsb_assembly_symmetry_axis = toTable(Schema.rcsb_assembly_symmetry_axis, model.sourceData.frame.categories.rcsb_assembly_symmetry_axis)
    } else {
        rcsb_assembly_symmetry_axis = toTable(Schema.rcsb_assembly_symmetry_axis, CifCategory.empty as any)
    }

    return _Database.ofTables('rcsb_assembly_symmetry', Schema, {
        rcsb_assembly_symmetry,
        rcsb_assembly_symmetry_cluster,
        rcsb_assembly_symmetry_cluster_member,
        rcsb_assembly_symmetry_axis
    })
}

const _Descriptor: ModelPropertyDescriptor = {
    isStatic: true,
    name: 'rcsb_assembly_symmetry',
    cifExport: {
        prefix: 'rcsb',
        categories: [
            PropertyWrapper.defaultInfoCategory<CifExportContext>('rcsb_assembly_symmetry_info', ctx => PropertyWrapper.createInfo()),
            getCategory('rcsb_assembly_symmetry'),
            getCategory('rcsb_assembly_symmetry_cluster'),
            getCategory('rcsb_assembly_symmetry_cluster_member'),
            getCategory('rcsb_assembly_symmetry_axis')
        ]
    }
}

export interface AssemblySymmetry {
    '@type': 'rcsb_assembly_symmetry',
    db: AssemblySymmetry.Database
    getSymmetries(assemblyIds: string[]): Table<AssemblySymmetry.Schema['rcsb_assembly_symmetry']>
    getClusters(symmetryId: number): Table<AssemblySymmetry.Schema['rcsb_assembly_symmetry_cluster']>
    getClusterMembers(clusterId: number): Table<AssemblySymmetry.Schema['rcsb_assembly_symmetry_cluster_member']>
    getAxes(symmetryId: number): Table<AssemblySymmetry.Schema['rcsb_assembly_symmetry_axis']>
}

export function AssemblySymmetry(db: AssemblySymmetry.Database): AssemblySymmetry {
    const f = db.rcsb_assembly_symmetry
    const c = db.rcsb_assembly_symmetry_cluster
    const cm = db.rcsb_assembly_symmetry_cluster_member
    const a = db.rcsb_assembly_symmetry_axis

    return {
        '@type': 'rcsb_assembly_symmetry',
        db,
        getSymmetries: (assemblyIds: string[]) => Table.pick(f, f._schema, i => assemblyIds.includes(f.assembly_id.value(i))),
        getClusters: (symmetryId: number) => Table.pick(c, c._schema, i => c.symmetry_id.value(i) === symmetryId),
        getClusterMembers: (clusterId: number) => Table.pick(cm, cm._schema, i => cm.cluster_id.value(i) === clusterId),
        getAxes: (symmetryId: number) => Table.pick(a, a._schema, i => a.symmetry_id.value(i) === symmetryId),
    }
}

const Client = new GraphQLClient(AssemblySymmetry.GraphQLEndpointURL, (url: string, type: 'string' | 'binary', body?: string) => ajaxGet({ url, type, body }) )

export namespace AssemblySymmetry {
    export function is(x: any): x is AssemblySymmetry {
        return x['@type'] === 'rcsb_assembly_symmetry'
    }
    export const GraphQLEndpointURL = 'http://rest-experimental.rcsb.org/graphql'
    export const Schema = {
        rcsb_assembly_symmetry_info: {
            updated_datetime_utc: Column.Schema.str
        },
        rcsb_assembly_symmetry: {
            /** Uniquely identifies a record in `rcsb_assembly_symmetry` */
            id: int,
            /**
             * A pointer to `pdbx_struct_assembly.id`.
             * The value 'deposited' refers to the coordinates as given in the file.
             * */
            assembly_id: str,
            /** Name and version of software used to calculate assembly symmetry */
            provenance: str,
            /** Type of protein symmetry */
            kind: Aliased<'GLOBAL' | 'LOCAL' | 'PSEUDO'>(str),
            /** Quantitative description of every individual subunit in a given protein */
            stoichiometry: List(',', x => x),
            /**
             * Symmetry symbol refers to point group or helical symmetry of identical subunits.
             * Contains point group symbol (e.g., C2, C5, D2, T, O, I) or H for helical symmetry.
             */
            symbol: str,
            /** Point group or helical symmetry */
            type: Aliased<'ASYMMETRIC' | 'CYCLIC' | 'DIHEDRAL' | 'HELICAL' | 'ICOSAHEDRAL' | 'OCTAHEDRAL' | 'TETRAHEDRAL'>(str),
            /**
             * Oligomeric state refers to a composition of subunits in quaternary structure.
             * Quaternary structure may be composed either exclusively of several copies of identical
             * subunits, in which case they are termed homo-oligomers, or alternatively by at least
             * one copy of different subunits (hetero-oligomers). Quaternary structure composed of
             * a single subunit is demoted as 'Monomer'.
             */
            oligomeric_state: str,
        },
        rcsb_assembly_symmetry_cluster: {
            /** Uniquely identifies a record in `rcsb_assembly_symmetry_cluster` */
            id: int,
            /** A pointer to `rcsb_assembly_symmetry.id` */
            symmetry_id: int,
            /** Average RMSD between members of a given cluster */
            avg_rmsd: float
        },
        rcsb_assembly_symmetry_cluster_member: {
            /** A pointer to `rcsb_assembly_symmetry_cluster.id` */
            cluster_id: int,
            /** The `label_asym_id` value of the member */
            asym_id: str,
            /** List of `pdbx_struct_oper_list_id` values of the member */
            pdbx_struct_oper_list_ids: List(',', x => x)
        },
        rcsb_assembly_symmetry_axis: {
            /** A pointer to `rcsb_assembly_symmetry.id` */
            symmetry_id: int,
            /**
             * The number of times (order of rotation) that a subunit can be repeated by a rotation
             * operation, being transformed into a new state indistinguishable from its starting state.
             */
            order: int,
            /** The x,y,z coordinate of the start point of a symmetry axis */
            start: Vector(3),
            /** The x,y,z coordinate of the end point of a symmetry axis */
            end: Vector(3)
        }
    }
    export type Schema = typeof Schema
    export interface Database extends _Database<Schema> {}

    export const Descriptor = _Descriptor;

    export async function attachFromCifOrAPI(model: Model, client: GraphQLClient = Client, ctx?: RuntimeContext) {
        if (model.customProperties.has(Descriptor)) return true;

        let db: Database
        let info = PropertyWrapper.tryGetInfoFromCif('rcsb_assembly_symmetry_info', model);
        if (info) {
            db = createDatabaseFromCif(model)
        } else {
            let result: AssemblySymmetryGraphQL.Query
            const variables: AssemblySymmetryGraphQL.Variables = { pdbId: model.label.toLowerCase() };
            try {
                result = await client.request<AssemblySymmetryGraphQL.Query>(ctx || RuntimeContext.Synchronous, query, variables);
            } catch (e) {
                console.error(e)
                return false;
            }
            if (!result || !result.assemblies) return false;

            db = createDatabaseFromJson(result.assemblies as ReadonlyArray<AssemblySymmetryGraphQL.Assemblies>)
        }

        model.customProperties.add(Descriptor);
        model._staticPropertyData.__RCSBAssemblySymmetry__ = AssemblySymmetry(db);
        return true;
    }

    export function createAttachTask(fetch: (url: string, type: 'string' | 'binary') => Task<string | Uint8Array>) {
        return (model: Model) => Task.create('RCSB Assembly Symmetry', async ctx => {
            if (get(model)) return true;

            return await attachFromCifOrAPI(model, new GraphQLClient(AssemblySymmetry.GraphQLEndpointURL, fetch), ctx)
        });
    }

    export function get(model: Model): AssemblySymmetry | undefined {
        return model._staticPropertyData.__RCSBAssemblySymmetry__;
    }
}