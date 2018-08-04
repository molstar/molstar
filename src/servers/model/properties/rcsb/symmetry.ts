/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GraphQLClient } from 'graphql-request'

import { RcsbSymmetry } from './graphql/types';
import query from './graphql/symmetry.gql';

import { Model, ModelPropertyDescriptor } from 'mol-model/structure';
import { CifWriter } from 'mol-io/writer/cif';
import { Database as _Database, Column, Table } from 'mol-data/db'
import { Category } from 'mol-io/writer/cif/encoder';
import { Tensor } from 'mol-math/linear-algebra';
import { CifExportContext } from 'mol-model/structure/export/mmcif';

const { str, int, float, Aliased, Vector, List } = Column.Schema;

function getInstance(name: keyof SymmetryAnnotation.Schema): (ctx: CifExportContext) => CifWriter.Category.Instance<any, any> {
    return function(ctx: CifExportContext) {
        const db = SymmetryAnnotation.get(ctx.model);
        return db ? Category.ofTable(db[name]) : CifWriter.Category.Empty;
    }
}

function getCategory(name: keyof SymmetryAnnotation.Schema) {
    return { name, instance: getInstance(name) }
}

function createDatabase(assemblies: ReadonlyArray<RcsbSymmetry.Assemblies>): SymmetryAnnotation.Database {
    const Schema = SymmetryAnnotation.Schema

    const featureRows: Table.Row<typeof Schema.symmetry_annotation_feature>[] = []
    const clusterRows: Table.Row<typeof Schema.symmetry_annotation_cluster>[] = []
    const axisRows: Table.Row<typeof Schema.symmetry_annotation_axis>[] = []

    let id = 0
    for (let i = 0, il = assemblies.length; i < il; ++i) {
        const a = assemblies[i]
        const assembly_id = (a.assembly_id!).toString()
        const source = a.rcsb_annotation_symmetry!.source!
        const symmetry_features = a.rcsb_annotation_symmetry!.symmetry_features!
        for (let j = 0, jl = symmetry_features.length; j < jl; ++j) {
            const f = symmetry_features[j]!
            featureRows.push({
                id,
                assembly_id,
                source,
                type: f.type!,
                stoichiometry_value: (f.stoichiometry!.value!) as string[],
                stoichiometry_description: f.stoichiometry!.description!
            })

            const clusters = f.clusters
            if (clusters) {
                for (let k = 0, kl = clusters.length; k < kl; ++k) {
                    const c = clusters[k]!
                    clusterRows.push({
                        feature_id: id,
                        avg_rmsd: c.avg_rmsd!,
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
                        order: a.order!,
                        start: a.start as Tensor.Data,
                        end: a.end as Tensor.Data
                    })
                }
            }

            id += 1
        }
    }

    return _Database.ofTables('symmetry_annotation', Schema, {
        symmetry_annotation_feature: Table.ofRows(Schema.symmetry_annotation_feature, featureRows),
        symmetry_annotation_cluster: Table.ofRows(Schema.symmetry_annotation_cluster, clusterRows),
        symmetry_annotation_axis: Table.ofRows(Schema.symmetry_annotation_axis, axisRows)
    })
}

const _Descriptor: ModelPropertyDescriptor = {
    isStatic: true,
    name: 'symmetry_annotation',
    cifExport: {
        categories: [
            getCategory('symmetry_annotation_feature'),
            getCategory('symmetry_annotation_cluster'),
            getCategory('symmetry_annotation_axis')
        ]
    }
}

const client = new GraphQLClient('http://rest-experimental.rcsb.org/graphql')

export namespace SymmetryAnnotation {
    export const Schema = {
        symmetry_annotation_feature: {
            id: int,
            assembly_id: str,
            source: str,
            type: Aliased<'GLOBAL' | 'LOCAL' | 'PSEUDO'>(str),
            stoichiometry_value: List(',', x => x),
            stoichiometry_description: str
        },
        symmetry_annotation_cluster: {
            feature_id: int,
            avg_rmsd: float,
            members: List(',', x => x)
        },
        symmetry_annotation_axis: {
            feature_id: int,
            order: int,
            start: Vector(3),
            end: Vector(3)
        }
    }
    export type Schema = typeof Schema
    export interface Database extends _Database<Schema> {}

    export const Descriptor = _Descriptor;

    export async function attachFromRCSB(model: Model) {
        if (model.customProperties.has(Descriptor)) return true;

        const variables: RcsbSymmetry.Variables = { pdbId: model.label.toLowerCase() };
        const result = await client.request<RcsbSymmetry.Query>(query, variables);
        if (!result || !result.assemblies) return false;

        const db: Database = createDatabase(result.assemblies as ReadonlyArray<RcsbSymmetry.Assemblies>)
        model.customProperties.add(Descriptor);
        model._staticPropertyData.__SymmetryAnnotation__ = db;

        return true;
    }

    export function get(model: Model): Database | undefined {
        return model._staticPropertyData.__SymmetryAnnotation__;
    }
}