/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ThemeDataContext } from 'mol-theme/theme';
import { ColorTheme, LocationColor } from 'mol-theme/color';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { Table } from 'mol-data/db';
import { AssemblySymmetry } from '../assembly-symmetry';
import { ColorScale, Color } from 'mol-util/color';
import { Unit, StructureElement, StructureProperties } from 'mol-model/structure';
import { Location } from 'mol-model/location';
import { ColorListName, ColorListOptions } from 'mol-util/color/scale';

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

function clusterMemberKey (asym_id: string, oper_list_ids: string[]) {
    return `${asym_id}-${oper_list_ids.join('x')}`
}

export const AssemblySymmetryClusterColorThemeParams = {
    list: PD.Select<ColorListName>('Viridis', ColorListOptions),
    symmetryId: PD.Select<number>(-1, []),
}
export type AssemblySymmetryClusterColorThemeParams = typeof AssemblySymmetryClusterColorThemeParams
export function getAssemblySymmetryClusterColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(AssemblySymmetryClusterColorThemeParams)

    if (ctx.structure && ctx.structure.models[0].customProperties.has(AssemblySymmetry.Descriptor)) {
        const assemblySymmetry = AssemblySymmetry.get(ctx.structure.models[0])!
        const assemblyName = ctx.structure.assemblyName
        const s = assemblySymmetry.db.rcsb_assembly_symmetry
        if (s._rowCount) {
            params.symmetryId.options = []
            for (let i = 0, il = s._rowCount; i < il; ++i) {
                if (s.assembly_id.value(i) === assemblyName) {
                    params.symmetryId.options.push([
                        s.id.value(i), `${s.symbol.value(i)} ${s.kind.value(i)}`
                    ])
                }
            }
            params.symmetryId.defaultValue = params.symmetryId.options[0][0]
        }
    }

    return params
}

export function AssemblySymmetryClusterColorTheme(ctx: ThemeDataContext, props: PD.Values<AssemblySymmetryClusterColorThemeParams>): ColorTheme<AssemblySymmetryClusterColorThemeParams> {
    let color: LocationColor = () => DefaultColor

    const { symmetryId } = props

    if (ctx.structure && ctx.structure.models[0].customProperties.has(AssemblySymmetry.Descriptor)) {
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
                        const asym_id = clusterMembers.asym_id.value(j)
                        const oper_list_ids = clusterMembers.pdbx_struct_oper_list_ids.value(j)
                        if (oper_list_ids.length === 0) oper_list_ids.push('1') // TODO hack assuming '1' is the id of the identity operator
                        clusterByMember.set(clusterMemberKey(asym_id, oper_list_ids), i)
                    }
                }

                const scale = ColorScale.create({ listOrName: props.list, domain: [ 0, clusters._rowCount - 1 ] })

                color = (location: Location): Color => {
                    if (StructureElement.isLocation(location)) {
                        const asym_id = getAsymId(location.unit)
                        const ns = location.unit.conformation.operator.name.split('-')
                        const oper_list_ids = ns.length === 2 ? ns[1].split('x') : []
                        const cluster = clusterByMember.get(clusterMemberKey(asym_id(location), oper_list_ids))
                        return cluster !== undefined ? scale.color(cluster) : DefaultColor
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
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && ctx.structure.models[0].customProperties.has(AssemblySymmetry.Descriptor)
}