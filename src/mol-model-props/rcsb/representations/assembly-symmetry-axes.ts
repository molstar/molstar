/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { VisualContext } from '../../../mol-repr/visual';
import { ThemeRegistryContext, Theme } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { StructureRepresentationProvider, StructureRepresentation, ComplexRepresentation, ComplexVisual } from '../../../mol-repr/structure/representation';
import { AssemblySymmetry } from '../assembly-symmetry';
import { Table } from '../../../mol-data/db';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Tensor } from '../../../mol-math/linear-algebra';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { addCylinder } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { VisualUpdateState } from '../../../mol-repr/util';
import { ComplexMeshVisual, ComplexMeshParams } from '../../../mol-repr/structure/complex-visual';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { EmptyLoci, createDataLoci, Loci, isDataLoci } from '../../../mol-model/loci';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { NullLocation } from '../../../mol-model/location';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { OrderedSet, Interval } from '../../../mol-data/int';
import { getSymmetrySelectParam, getAssemblyIds } from '../util';

export const AssemblySymmetryAxesParams = {
    symmetryId: getSymmetrySelectParam(),
    sizeFactor: PD.Numeric(0.4, { min: 0, max: 3, step: 0.01 }),

    ...ComplexMeshParams,
    radialSegments: PD.Numeric(16, { min: 3, max: 56, step: 1 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }),
}
export type AssemblySymmetryAxesParams = typeof AssemblySymmetryAxesParams
export function getAssemblySymmetryAxesParams(ctx: ThemeRegistryContext, structure: Structure) {
    const params = PD.clone(AssemblySymmetryAxesParams)
    params.symmetryId = getSymmetrySelectParam(structure)
    params.unitKinds.isHidden = true
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

export function AssemblySymmetryAxesVisual(materialId: number): ComplexVisual<AssemblySymmetryAxesParams> {
    return ComplexMeshVisual<AssemblySymmetryAxesParams>({
        defaultProps: PD.getDefaultValues(AssemblySymmetryAxesParams),
        createGeometry: createAssemblySymmetryAxesMesh,
        createLocationIterator,
        getLoci,
        eachLocation: eachAxisLocation,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<AssemblySymmetryAxesParams>, currentProps: PD.Values<AssemblySymmetryAxesParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail ||
                newProps.symmetryId !== currentProps.symmetryId
            )
        }
    }, materialId)
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

function eachAxisLocation(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
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

    // check if structure.units operators have symmetry.assembly_id
    if (!getAssemblyIds(structure.units).includes(symmetry.assembly_id)) return Mesh.createEmpty(mesh)

    const axes = assemblySymmetry.db.rcsb_assembly_symmetry_axis
    const vectorSpace = AssemblySymmetry.Schema.rcsb_assembly_symmetry_axis.start.space;

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
    }
    return MeshBuilder.getMesh(builderState)
}