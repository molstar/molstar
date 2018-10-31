/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssemblySymmetry } from 'mol-model-props/rcsb/symmetry';
import { Table } from 'mol-data/db';
import { Color, ColorScale } from 'mol-util/color';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { Tensor } from 'mol-math/linear-algebra';
import { addSphere } from 'mol-geo/geometry/mesh/builder/sphere';
import { addCylinder } from 'mol-geo/geometry/mesh/builder/cylinder';
import { Shape } from 'mol-model/shape';
import { ColorTheme } from 'mol-theme/color';
import { Location } from 'mol-model/location';
import { StructureElement, Unit, StructureProperties } from 'mol-model/structure';

export function getAxesShape(symmetryId: number, assemblySymmetry: AssemblySymmetry) {
    const s = assemblySymmetry.db.rcsb_assembly_symmetry
    const symmetry = Table.pickRow(s, i => s.id.value(i) === symmetryId)
    if (!symmetry) return

    const axes = assemblySymmetry.getAxes(symmetryId)
    if (!axes._rowCount) return

    const vectorSpace = AssemblySymmetry.Schema.rcsb_assembly_symmetry_axis.start.space;

    const colors: Color[] = []
    const labels: string[] = []

    const radius = 0.4
    const cylinderProps = { radiusTop: radius, radiusBottom: radius }
    const meshBuilder = MeshBuilder.create(256, 128)

    for (let i = 0, il = axes._rowCount; i < il; ++i) {
        const start = Tensor.toVec3(vectorSpace, axes.start.value(i))
        const end = Tensor.toVec3(vectorSpace, axes.end.value(i))
        meshBuilder.setGroup(i)
        addSphere(meshBuilder, start, radius, 2)
        addSphere(meshBuilder, end, radius, 2)
        addCylinder(meshBuilder, start, end, 1, cylinderProps)
        colors.push(Color(0xCCEE11))
        labels.push(`Axis ${i + 1} for ${symmetry.kind} ${symmetry.type.toLowerCase()} symmetry`)
    }
    const mesh = meshBuilder.getMesh()
    const shape = Shape.create('Axes', mesh, colors, labels)
    return shape
}

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

export function getClusterColorTheme(symmetryId: number, assemblySymmetry: AssemblySymmetry): ColorTheme {
    const DefaultColor = Color(0xCCCCCC)
    const s = assemblySymmetry.db.rcsb_assembly_symmetry
    const symmetry = Table.pickRow(s, i => s.id.value(i) === symmetryId)
    if (!symmetry) return { granularity: 'uniform', color: () => DefaultColor }

    const clusters = assemblySymmetry.getClusters(symmetryId)
    if (!clusters._rowCount) return { granularity: 'uniform', color: () => DefaultColor }

    const clusterByMember = new Map<string, number>()
    for (let i = 0, il = clusters._rowCount; i < il; ++i) {
        const clusterMembers = assemblySymmetry.getClusterMembers(clusters.id.value(i))
        for (let j = 0, jl = clusterMembers._rowCount; j < jl; ++j) {
            const asym_id = clusterMembers.asym_id.value(j)
            const oper_list_ids = clusterMembers.pdbx_struct_oper_list_ids.value(j)
            clusterByMember.set(clusterMemberKey(asym_id, oper_list_ids), i)
        }
    }
    const scale = ColorScale.create({ domain: [ 0, clusters._rowCount - 1 ] })

    return {
        granularity: 'instance',
        color: (location: Location): Color => {
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