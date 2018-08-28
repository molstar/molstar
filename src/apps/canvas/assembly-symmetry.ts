/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssemblySymmetry } from 'mol-model-props/rcsb/symmetry';
import { Table } from 'mol-data/db';
import { Color } from 'mol-util/color';
import { MeshBuilder } from 'mol-geo/mesh/mesh-builder';
import { Tensor } from 'mol-math/linear-algebra';
import { addSphere } from 'mol-geo/mesh/builder/sphere';
import { addCylinder } from 'mol-geo/mesh/builder/cylinder';
import { Shape } from 'mol-model/shape';

export function getAxesShape(featureId: number, assemblySymmetry: AssemblySymmetry) {
    const f = assemblySymmetry.db.rcsb_assembly_symmetry_feature
    const feature = Table.pickRow(f, i => f.id.value(i) === featureId)
    if (!feature) return

    const axes = assemblySymmetry.getAxes(featureId)
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
        labels.push(`Axis ${i + 1} for ${feature.symmetry_value} ${feature.type.toLowerCase()} symmetry`)
    }
    const mesh = meshBuilder.getMesh()
    const shape = Shape.create('Axes', mesh, colors, labels)
    return shape
}

export function getClusterColorTheme(featureId: number, assemblySymmetry: AssemblySymmetry) {
    const f = assemblySymmetry.db.rcsb_assembly_symmetry_feature
    const feature = Table.pickRow(f, i => f.id.value(i) === featureId)
    if (!feature) return

    const clusters = assemblySymmetry.getClusters(featureId)
    if (!clusters._rowCount) return

    for (let i = 0, il = clusters._rowCount; i < il; ++i) {
        console.log(clusters.members.value(i), clusters.avg_rmsd.value(i), feature.stoichiometry_value, feature.stoichiometry_description)
    }
}