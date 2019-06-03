/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ThemeDataContext } from '../../../mol-theme/theme';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition'
import { Table } from '../../../mol-data/db';
import { AssemblySymmetry } from '../assembly-symmetry';
import { ColorScale, Color } from '../../../mol-util/color';
import { Unit, StructureElement, StructureProperties } from '../../../mol-model/structure';
import { Location } from '../../../mol-model/location';
import { ColorListName, ColorListOptions } from '../../../mol-util/color/scale';
import { getSymmetrySelectParam } from '../util';

const DefaultColor = Color(0xCCCCCC)

function getAsymId(unit: Unit): StructureElement.Property<string> {
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            return StructureProperties.chain.label_asym_id
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return StructureProperties.coarse.asym_id
    }
}

function clusterMemberKey(assemblyId: string, asymId: string, operList: string[]) {
    return `${assemblyId}-${asymId}-${operList.join('|')}`
}

export const AssemblySymmetryClusterColorThemeParams = {
    list: PD.Select<ColorListName>('Viridis', ColorListOptions),
    symmetryId: getSymmetrySelectParam(),
}
export type AssemblySymmetryClusterColorThemeParams = typeof AssemblySymmetryClusterColorThemeParams
export function getAssemblySymmetryClusterColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(AssemblySymmetryClusterColorThemeParams)
    params.symmetryId = getSymmetrySelectParam(ctx.structure)
    return params
}

export function AssemblySymmetryClusterColorTheme(ctx: ThemeDataContext, props: PD.Values<AssemblySymmetryClusterColorThemeParams>): ColorTheme<AssemblySymmetryClusterColorThemeParams> {
    let color: LocationColor = () => DefaultColor

    const { symmetryId } = props

    if (ctx.structure && !ctx.structure.isEmpty && ctx.structure.models[0].customProperties.has(AssemblySymmetry.Descriptor)) {
        const assemblySymmetry = AssemblySymmetry.get(ctx.structure.models[0])!

        const s = assemblySymmetry.db.rcsb_assembly_symmetry
        const symmetry = Table.pickRow(s, i => s.id.value(i) === symmetryId)
        if (symmetry) {

            const clusters = assemblySymmetry.getClusters(symmetryId)
            if (clusters._rowCount) {

                const clusterByMember = new Map<string, number>()
                for (let i = 0, il = clusters._rowCount; i < il; ++i) {
                    const clusterMembers = assemblySymmetry.getClusterMembers(clusters.id.value(i))
                    for (let j = 0, jl = clusterMembers._rowCount; j < jl; ++j) {
                        const asymId = clusterMembers.asym_id.value(j)
                        const operList = clusterMembers.pdbx_struct_oper_list_ids.value(j)
                        if (operList.length === 0) operList.push('1') // TODO hack assuming '1' is the id of the identity operator
                        clusterByMember.set(clusterMemberKey(symmetry.assembly_id, asymId, operList), i)
                    }
                }

                const scale = ColorScale.create({ listOrName: props.list, domain: [ 0, clusters._rowCount - 1 ] })

                color = (location: Location): Color => {
                    if (StructureElement.isLocation(location)) {
                        const { assembly } = location.unit.conformation.operator
                        if (assembly && assembly.id === symmetry.assembly_id) {
                            const asymId = getAsymId(location.unit)(location)
                            const cluster = clusterByMember.get(clusterMemberKey(assembly.id, asymId, assembly.operList))
                            return cluster !== undefined ? scale.color(cluster) : DefaultColor
                        }
                    }
                    return DefaultColor
                }
            }
        }
    }

    return {
        factory: AssemblySymmetryClusterColorTheme,
        granularity: 'instance',
        color,
        props,
        description: 'Assigns chain colors according to assembly symmetry cluster membership.',
    }
}

export const AssemblySymmetryClusterColorThemeProvider: ColorTheme.Provider<AssemblySymmetryClusterColorThemeParams> = {
    label: 'RCSB Assembly Symmetry Cluster',
    factory: AssemblySymmetryClusterColorTheme,
    getParams: getAssemblySymmetryClusterColorThemeParams,
    defaultValues: PD.getDefaultValues(AssemblySymmetryClusterColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && !ctx.structure.isEmpty && ctx.structure.models[0].customProperties.has(AssemblySymmetry.Descriptor)
}