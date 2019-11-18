/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, StructureElement } from '../../mol-model/structure';
import { ShapeRepresentation } from '../../mol-repr/shape/representation';
import { Shape } from '../../mol-model/shape';
import { ColorNames } from '../../mol-util/color/names';
import { RuntimeContext } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Vec3, Mat4 } from '../../mol-math/linear-algebra';
import Matrix from '../../mol-math/linear-algebra/matrix/matrix';
import { PrincipalAxes } from '../../mol-math/linear-algebra/matrix/principal-axes';
import { createCage } from '../../mol-geo/primitive/cage';
import { stringToWords } from '../../mol-util/string';
import { structureElementStatsLabel } from '../../mol-theme/label';

const tmpMatrixPos = Vec3.zero()
function getPositions(structure: Structure) {
    const positions = new Float32Array(structure.elementCount * 3)
    for (let i = 0, m = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i]
        const { elements } = unit
        const pos = unit.conformation.position
        for (let j = 0, jl = elements.length; j < jl; ++j) {
            pos(elements[j], tmpMatrixPos)
            Vec3.toArray(tmpMatrixPos, positions, m + j * 3)
        }
        m += elements.length * 3
    }
    return positions
}

interface OrientationData {
    label: string
    principalAxes: PrincipalAxes
    projectedScale: { d1a: number, d2a: number, d3a: number, d1b: number, d2b: number, d3b: number }
}

const OrientationVisuals = { 'principal-axes': '', 'oriented-box': '' }
type OrientationVisualName = keyof typeof OrientationVisuals
const OrientationVisualOptions = Object.keys(OrientationVisuals).map(name => [name, stringToWords(name)] as [OrientationVisualName, string])

export const OrientationParams = {
    ...Mesh.Params,
    visuals: PD.MultiSelect<OrientationVisualName>(['oriented-box'], OrientationVisualOptions),
    orientationColor: PD.Color(ColorNames.orange),
    orientationScale: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 })
}
export type OrientationParams = typeof OrientationParams
export type OrientationProps = PD.Values<OrientationParams>

enum VisualGroup {
    PrincipalAxes = 1,
    OrientedBox = 2
}

function getVolume(data: OrientationData) {
    const { d1a, d2a, d3a, d1b, d2b, d3b } = data.projectedScale
    return (d1a - d1b) * (d2a - d2b) * (d3a - d3b)
}

function buildPrincipalAxes(state: MeshBuilder.State, data: OrientationData, props: OrientationProps) {
    const vertices = new Float32Array(6 * 3)
    const edges = new Uint8Array([0, 1, 2, 3, 4, 5])
    Vec3.toArray(data.principalAxes.begA, vertices, 0)
    Vec3.toArray(data.principalAxes.endA, vertices, 3)
    Vec3.toArray(data.principalAxes.begB, vertices, 6)
    Vec3.toArray(data.principalAxes.endB, vertices, 9)
    Vec3.toArray(data.principalAxes.begC, vertices, 12)
    Vec3.toArray(data.principalAxes.endC, vertices, 15)

    const matrix = Mat4.fromTranslation(Mat4(), Vec3.inverse(Vec3(), data.principalAxes.center))

    const cage = createCage(vertices, edges)
    const radius = (Math.cbrt(getVolume(data)) / 300) * props.orientationScale
    state.currentGroup = VisualGroup.PrincipalAxes
    MeshBuilder.addCage(state, matrix, cage, radius, 2, 20)
}

const tmpBoxVec = Vec3()
function buildOrientedBox(state: MeshBuilder.State, data: OrientationData, props: OrientationProps) {
    const { center, normVecA, normVecB, normVecC } = data.principalAxes
    const { d1a, d2a, d3a, d1b, d2b, d3b } = data.projectedScale

    const vertices = new Float32Array(8 * 3)
    const edges = new Uint8Array([
        0, 1, 0, 3, 0, 6, 1, 2, 1, 7, 2, 3,
        2, 4, 3, 5, 4, 5, 4, 7, 5, 6, 6, 7
    ])

    let offset = 0
    const addCornerHelper = function (d1: number, d2: number, d3: number) {
        Vec3.copy(tmpBoxVec, center)
        Vec3.scaleAndAdd(tmpBoxVec, tmpBoxVec, normVecA, d1)
        Vec3.scaleAndAdd(tmpBoxVec, tmpBoxVec, normVecB, d2)
        Vec3.scaleAndAdd(tmpBoxVec, tmpBoxVec, normVecC, d3)
        Vec3.toArray(tmpBoxVec, vertices, offset)
        offset += 3
    }
    addCornerHelper(d1a, d2a, d3a)
    addCornerHelper(d1a, d2a, d3b)
    addCornerHelper(d1a, d2b, d3b)
    addCornerHelper(d1a, d2b, d3a)
    addCornerHelper(d1b, d2b, d3b)
    addCornerHelper(d1b, d2b, d3a)
    addCornerHelper(d1b, d2a, d3a)
    addCornerHelper(d1b, d2a, d3b)

    const matrix = Mat4.fromTranslation(Mat4(), Vec3.inverse(Vec3(), data.principalAxes.center))

    const cage = createCage(vertices, edges)
    const radius = (Math.cbrt(getVolume(data)) / 300) * props.orientationScale
    state.currentGroup = VisualGroup.OrientedBox
    MeshBuilder.addCage(state, matrix, cage, radius, 2, 20)
}

function getOrientationMesh(data: OrientationData, props: OrientationProps, mesh?: Mesh) {
    const state = MeshBuilder.createState(256, 128, mesh)

    if (props.visuals.includes('principal-axes')) buildPrincipalAxes(state, data, props)
    if (props.visuals.includes('oriented-box')) buildOrientedBox(state, data, props)

    return MeshBuilder.getMesh(state)
}

function getLabel(structure: Structure) {
    const loci = Structure.toStructureElementLoci(structure)
    const remappedLoci = StructureElement.Loci.remap(loci, structure.root)
    return structureElementStatsLabel(StructureElement.Stats.ofLoci(remappedLoci), true)
}

export async function getStructureOrientationRepresentation(ctx: RuntimeContext, structure: Structure, params: OrientationProps, prev?: ShapeRepresentation<OrientationData, Mesh, Mesh.Params>) {
    const repr = prev || ShapeRepresentation(getOrientationShape, Mesh.Utils);
    const label = getLabel(structure)
    const positions = getPositions(structure)
    const principalAxes = PrincipalAxes.ofPoints(Matrix.fromArray(positions, 3, structure.elementCount))
    const projectedScale = PrincipalAxes.getProjectedScale(positions, principalAxes)
    const data: OrientationData = { label, principalAxes, projectedScale }
    await repr.createOrUpdate(params, data).runInContext(ctx);
    return repr;
}

function getOrientationLabel(data: OrientationData, groupId: number): string {
    switch (groupId) {
        case VisualGroup.PrincipalAxes: return 'Principal Axes'
        case VisualGroup.OrientedBox: return 'Oriented Box'
    }
    return 'Unknown Orientation Visual'
}

function getOrientationShape(ctx: RuntimeContext, data: OrientationData, props: OrientationProps, shape?: Shape<Mesh>) {
    const geo = getOrientationMesh(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return `${getOrientationLabel(data, groupId)} of ${data.label}`
    }
    return Shape.create('Principal Axes', data, geo, () => props.orientationColor, () => 1, getLabel)
}