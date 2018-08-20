/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GraphQLClient } from 'graphql-request'

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

const { str, int, float, Aliased, Vector, List } = Column.Schema;

function getInstance(name: keyof AssemblySymmetry.Schema): (ctx: CifExportContext) => CifWriter.Category.Instance<any, any> {
    return function(ctx: CifExportContext) {
        const db = AssemblySymmetry.get(ctx.model);
        return db ? Category.ofTable(db[name]) : CifWriter.Category.Empty;
    }
}

function getCategory(name: keyof AssemblySymmetry.Schema) {
    return { name, instance: getInstance(name) }
}

function createDatabase(assemblies: ReadonlyArray<AssemblySymmetryGraphQL.Assemblies>): AssemblySymmetry.Database {
    const Schema = AssemblySymmetry.Schema

    const featureRows: Table.Row<typeof Schema.rcsb_assembly_symmetry_feature>[] = []
    const clusterRows: Table.Row<typeof Schema.rcsb_assembly_symmetry_cluster>[] = []
    const axisRows: Table.Row<typeof Schema.rcsb_assembly_symmetry_axis>[] = []

    let id = 0
    for (let i = 0, il = assemblies.length; i < il; ++i) {
        const { assembly_id: _assembly_id, rcsb_assembly_symmetry } = assemblies[i]
        if (!rcsb_assembly_symmetry) continue
        const assembly_id = _assembly_id.toString() // TODO type will be fixed to string upstream
        const source = rcsb_assembly_symmetry.source
        const symmetry_features = rcsb_assembly_symmetry.symmetry_features
        for (let j = 0, jl = symmetry_features.length; j < jl; ++j) {
            const f = symmetry_features[j]! // TODO upstream, array members should not be nullable
            featureRows.push({
                id,
                assembly_id,
                source,
                type: f.type,
                stoichiometry_value: f.stoichiometry.value as string[],  // TODO upstream, array members should not be nullable
                stoichiometry_description: f.stoichiometry.description,
                symmetry_value: f.symmetry.value,
                symmetry_description: f.symmetry.description
            })

            const clusters = f.clusters
            if (clusters) {
                for (let k = 0, kl = clusters.length; k < kl; ++k) {
                    const c = clusters[k]! // TODO upstream, array members should not be nullable
                    clusterRows.push({
                        feature_id: id,
                        avg_rmsd: c.avg_rmsd || 0, // TODO upstream, should not be nullable, or???
                        members: c.members as string[]
                    })
                }
            }

            const axes = f.symmetry_axes
            if (axes) {
                for (let k = 0, kl = axes.length; k < kl; ++k) {
                    const a = axes[k]!
                    axisRows.push({
                        feature_id: id,
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
        rcsb_assembly_symmetry_feature: Table.ofRows(Schema.rcsb_assembly_symmetry_feature, featureRows),
        rcsb_assembly_symmetry_cluster: Table.ofRows(Schema.rcsb_assembly_symmetry_cluster, clusterRows),
        rcsb_assembly_symmetry_axis: Table.ofRows(Schema.rcsb_assembly_symmetry_axis, axisRows)
    })
}

const _Descriptor: ModelPropertyDescriptor = {
    isStatic: true,
    name: 'assembly_symmetry',
    cifExport: {
        prefix: 'rcsb',
        categories: [
            getCategory('rcsb_assembly_symmetry_feature'),
            getCategory('rcsb_assembly_symmetry_cluster'),
            getCategory('rcsb_assembly_symmetry_axis')
        ]
    }
}

const client = new GraphQLClient('http://rest-experimental.rcsb.org/graphql')

export namespace AssemblySymmetry {
    export const Schema = {
        rcsb_assembly_symmetry_feature: {
            /** Uniquely identifies a record in `rcsb_assembly_symmetry_feature` */
            id: int,
            /** A pointer to `pdbx_struct_assembly.id` */
            assembly_id: str,
            /** Name and version of software used to calculate protein symmetry */
            source: str,
            /** Type of protein symmetry */
            type: Aliased<'GLOBAL' | 'LOCAL' | 'PSEUDO'>(str),
            /** Quantitative description of every individual subunit in a given protein */
            stoichiometry_value: List(',', x => x),
            /** Oligomeric state for a given composition of subunits */
            stoichiometry_description: str,
            /** Point group symmetry value */
            symmetry_value: str,
            /** Point group symmetry description */
            symmetry_description: str
        },
        rcsb_assembly_symmetry_cluster: {
            /** A pointer to `rcsb_assembly_symmetry_feature.id` */
            feature_id: int,
            /** Average RMSD between members of a given cluster */
            avg_rmsd: float,
            /** List of `auth_label_id` values  */
            members: List(',', x => x)
        },
        rcsb_assembly_symmetry_axis: {
            /** A pointer to `rcsb_assembly_symmetry_feature.id` */
            feature_id: int,
            /** The order of the symmetry axis */
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

    export async function attachFromCifOrAPI(model: Model) {
        if (model.customProperties.has(Descriptor)) return true;

        let db: Database

        if (model.sourceData.kind === 'mmCIF' && model.sourceData.frame.categoryNames.includes('rcsb_assembly_symmetry_feature')) {
            const rcsb_assembly_symmetry_feature = toTable(Schema.rcsb_assembly_symmetry_feature, model.sourceData.frame.categories.rcsb_assembly_symmetry_feature)

            let rcsb_assembly_symmetry_cluster
            if (model.sourceData.frame.categoryNames.includes('rcsb_assembly_symmetry_cluster')) {
                rcsb_assembly_symmetry_cluster = toTable(Schema.rcsb_assembly_symmetry_cluster, model.sourceData.frame.categories.rcsb_assembly_symmetry_cluster)
            } else {
                rcsb_assembly_symmetry_cluster = CifCategory.empty
            }

            let rcsb_assembly_symmetry_axis
            if (model.sourceData.frame.categoryNames.includes('rcsb_assembly_symmetry_axis')) {
                rcsb_assembly_symmetry_axis = toTable(Schema.rcsb_assembly_symmetry_axis, model.sourceData.frame.categories.rcsb_assembly_symmetry_axis)
            } else {
                rcsb_assembly_symmetry_axis = CifCategory.empty
            }

            db = _Database.ofTables('assembly_symmetry', Schema, {
                rcsb_assembly_symmetry_feature,
                rcsb_assembly_symmetry_cluster,
                rcsb_assembly_symmetry_axis
            })
        } else {
            let result: AssemblySymmetryGraphQL.Query
            const variables: AssemblySymmetryGraphQL.Variables = { pdbId: model.label.toLowerCase() };
            try {
                result = await client.request<AssemblySymmetryGraphQL.Query>(query, variables);
            } catch (e) {
                console.error(e)
                return false;
            }
            if (!result || !result.assemblies) return false;

            db = createDatabase(result.assemblies as ReadonlyArray<AssemblySymmetryGraphQL.Assemblies>)
        }

        model.customProperties.add(Descriptor);
        model._staticPropertyData.__AssemblySymmetry__ = db;
        return true;
    }

    export function get(model: Model): Database | undefined {
        return model._staticPropertyData.__AssemblySymmetry__;
    }
}