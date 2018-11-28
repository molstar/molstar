/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from 'mol-util/param-definition';
import { RepresentationParamsGetter, RepresentationContext, VisualContext } from 'mol-repr/representation';
import { ThemeRegistryContext, Theme } from 'mol-theme/theme';
import { Structure } from 'mol-model/structure';
import { StructureRepresentationProvider, StructureRepresentation, ComplexRepresentation, ComplexVisual } from 'mol-repr/structure/representation';
import { AssemblySymmetry } from '../assembly-symmetry';
import { Table } from 'mol-data/db';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { Tensor } from 'mol-math/linear-algebra';
import { addSphere } from 'mol-geo/geometry/mesh/builder/sphere';
import { addCylinder } from 'mol-geo/geometry/mesh/builder/cylinder';
import { VisualUpdateState } from 'mol-repr/util';
import { ComplexMeshVisual, ComplexMeshParams } from 'mol-repr/structure/complex-visual';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { EmptyLoci, createDataLoci, Loci, isDataLoci } from 'mol-model/loci';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { NullLocation } from 'mol-model/location';
import { PickingId } from 'mol-geo/geometry/picking';
import { OrderedSet, Interval } from 'mol-data/int';

export const AssemblySymmetryAxesParams = {
    ...ComplexMeshParams,
    sizeFactor: PD.Numeric(0.4, { min: 0, max: 3, step: 0.01 }),
    radialSegments: PD.Numeric(16, { min: 3, max: 56, step: 1 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }),
    symmetryId: PD.Select<number>(-1, []),
}
export type AssemblySymmetryAxesParams = typeof AssemblySymmetryAxesParams
export function getAssemblySymmetryAxesParams(ctx: ThemeRegistryContext, structure: Structure) {
    const params = PD.clone(AssemblySymmetryAxesParams)

    if (structure.models[0].customProperties.has(AssemblySymmetry.Descriptor)) {
        const assemblySymmetry = AssemblySymmetry.get(structure.models[0])!
        const assemblyName = structure.assemblyName
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

export type AssemblySymmetryAxesRepresentation = StructureRepresentation<AssemblySymmetryAxesParams>
export function AssemblySymmetryAxesRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, AssemblySymmetryAxesParams>): AssemblySymmetryAxesRepresentation {
    return ComplexRepresentation('RCSB Assembly Symmetry Axes', ctx, getParams, AssemblySymmetryAxesVisual)
}

export const AssemblySymmetryAxesRepresentationProvider: StructureRepresentationProvider<AssemblySymmetryAxesParams> = {
    label: 'RCSB Assembly Symmetry Axes',
    description: 'Displays assembly symmetry axes.',
    factory: AssemblySymmetryAxesRepresentation,
    getParams: getAssemblySymmetryAxesParams,
    defaultValues: PD.getDefaultValues(AssemblySymmetryAxesParams),
    defaultColorTheme: 'uniform',
    defaultSizeTheme: 'uniform'
}

//

export function AssemblySymmetryAxesVisual(): ComplexVisual<AssemblySymmetryAxesParams> {
    return ComplexMeshVisual<AssemblySymmetryAxesParams>({
        defaultProps: PD.getDefaultValues(AssemblySymmetryAxesParams),
        createGeometry: createAssemblySymmetryAxesMesh,
        createLocationIterator,
        getLoci,
        mark,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<AssemblySymmetryAxesParams>, currentProps: PD.Values<AssemblySymmetryAxesParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail ||
                newProps.symmetryId !== currentProps.symmetryId
            )
        }
    })
}

function createLocationIterator(structure: Structure) {
    const assemblySymmetry = AssemblySymmetry.get(structure.models[0])
    const groupCount = assemblySymmetry ? assemblySymmetry.db.rcsb_assembly_symmetry_axis._rowCount : 0
    return LocationIterator(groupCount, 1, () => NullLocation)
}

function getLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const assemblySymmetry = AssemblySymmetry.get(structure.models[0])
        if (assemblySymmetry) {
            return createDataLoci(assemblySymmetry, 'axes', OrderedSet.ofSingleton(groupId))
        }
    }
    return EmptyLoci
}

function mark(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false
    if (!isDataLoci(loci) || loci.tag !== 'axes') return false
    const assemblySymmetry = AssemblySymmetry.get(structure.models[0])
    if (!assemblySymmetry || loci.data !== assemblySymmetry) return false
    OrderedSet.forEach(loci.indices, v => { if (apply(Interval.ofSingleton(v))) changed = true })
    return changed
}

export function createAssemblySymmetryAxesMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<AssemblySymmetryAxesParams>, mesh?: Mesh) {

    const { symmetryId, sizeFactor } = props

    const assemblySymmetry = AssemblySymmetry.get(structure.models[0])
    if (!assemblySymmetry) return Mesh.createEmpty(mesh)

    const s = assemblySymmetry.db.rcsb_assembly_symmetry
    const symmetry = Table.pickRow(s, i => s.id.value(i) === symmetryId)
    if (!symmetry) return Mesh.createEmpty(mesh)

    // symmetry.assembly_id not available for structure.assemblyName
    if (symmetry.assembly_id !== structure.assemblyName) return Mesh.createEmpty(mesh)

    const axes = assemblySymmetry.db.rcsb_assembly_symmetry_axis
    const vectorSpace = AssemblySymmetry.Schema.rcsb_assembly_symmetry_axis.start.space;
    // const colors: Color[] = []
    // const labels: string[] = []

    const radius = 1 * sizeFactor
    const cylinderProps = { radiusTop: radius, radiusBottom: radius }
    const builderState = MeshBuilder.createState(256, 128, mesh)

     for (let i = 0, il = axes._rowCount; i < il; ++i) {
        if (axes.symmetry_id.value(i) !== symmetryId) continue

        const start = Tensor.toVec3(vectorSpace, axes.start.value(i))
        const end = Tensor.toVec3(vectorSpace, axes.end.value(i))
        builderState.currentGroup = i
        addSphere(builderState, start, radius, 2)
        addSphere(builderState, end, radius, 2)
        addCylinder(builderState, start, end, 1, cylinderProps)
        // colors.push(Color(0xCCEE11))
        // labels.push(`Axis ${i + 1} for ${symmetry.kind} ${symmetry.type.toLowerCase()} symmetry`)
    }
    return MeshBuilder.getMesh(builderState)
}