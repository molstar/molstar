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
import { PrincipalAxes } from '../../mol-math/linear-algebra/matrix/principal-axes';
import { createCage } from '../../mol-geo/primitive/cage';
import { stringToWords } from '../../mol-util/string';
import { structureElementStatsLabel } from '../../mol-theme/label';
import { Axes3D } from '../../mol-math/geometry';

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
}

const OrientationVisuals = { 'principal-axes': '', 'oriented-box': '' }
type OrientationVisualName = keyof typeof OrientationVisuals
const OrientationVisualOptions = Object.keys(OrientationVisuals).map(name => [name, stringToWords(name)] as [OrientationVisualName, string])

export const OrientationParams = {
    ...Mesh.Params,
    visuals: PD.MultiSelect<OrientationVisualName>(['oriented-box'], OrientationVisualOptions),
    color: PD.Color(ColorNames.orange),
    scale: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 })
}
export type OrientationParams = typeof OrientationParams
export type OrientationProps = PD.Values<OrientationParams>

enum VisualGroup {
    PrincipalAxes = 1,
    OrientedBox = 2
}

const tmpAxesVec = Vec3()
function buildMomentsAxes(state: MeshBuilder.State, data: OrientationData, props: OrientationProps) {
    const { origin, dirA, dirB, dirC } = data.principalAxes.momentsAxes

    const vertices = new Float32Array(6 * 3)
    const edges = new Uint8Array([0, 1, 2, 3, 4, 5])
    Vec3.add(tmpAxesVec, origin, dirA)
    Vec3.toArray(Vec3.add(tmpAxesVec, origin, dirA), vertices, 0)
    Vec3.toArray(Vec3.sub(tmpAxesVec, origin, dirA), vertices, 3)
    Vec3.toArray(Vec3.add(tmpAxesVec, origin, dirB), vertices, 6)
    Vec3.toArray(Vec3.sub(tmpAxesVec, origin, dirB), vertices, 9)
    Vec3.toArray(Vec3.add(tmpAxesVec, origin, dirC), vertices, 12)
    Vec3.toArray(Vec3.sub(tmpAxesVec, origin, dirC), vertices, 15)

    const matrix = Mat4.fromTranslation(Mat4(), Vec3.negate(Vec3(), origin))

    const cage = createCage(vertices, edges)
    const volume = Axes3D.volume(data.principalAxes.boxAxes)
    const radius = (Math.cbrt(volume) / 300) * props.scale
    state.currentGroup = VisualGroup.PrincipalAxes
    MeshBuilder.addCage(state, matrix, cage, radius, 2, 20)
}

const tmpBoxVecCorner = Vec3()
const tmpBoxVecA = Vec3()
const tmpBoxVecB = Vec3()
const tmpBoxVecC = Vec3()
function buildOrientedBox(state: MeshBuilder.State, data: OrientationData, props: OrientationProps) {
    const { origin, dirA, dirB, dirC } = data.principalAxes.boxAxes
    const negDirA = Vec3.negate(tmpBoxVecA, dirA)
    const negDirB = Vec3.negate(tmpBoxVecB, dirB)
    const negDirC = Vec3.negate(tmpBoxVecC, dirC)

    const vertices = new Float32Array(8 * 3)
    const edges = new Uint8Array([
        0, 1, 0, 3, 0, 6, 1, 2, 1, 7, 2, 3,
        2, 4, 3, 5, 4, 5, 4, 7, 5, 6, 6, 7
    ])

    let offset = 0
    const addCornerHelper = function (v1: Vec3, v2: Vec3, v3: Vec3) {
        Vec3.copy(tmpBoxVecCorner, origin)
        Vec3.add(tmpBoxVecCorner, tmpBoxVecCorner, v1)
        Vec3.add(tmpBoxVecCorner, tmpBoxVecCorner, v2)
        Vec3.add(tmpBoxVecCorner, tmpBoxVecCorner, v3)
        Vec3.toArray(tmpBoxVecCorner, vertices, offset)
        offset += 3
    }
    addCornerHelper(dirA, dirB, dirC)
    addCornerHelper(dirA, dirB, negDirC)
    addCornerHelper(dirA, negDirB, negDirC)
    addCornerHelper(dirA, negDirB, dirC)
    addCornerHelper(negDirA, negDirB, negDirC)
    addCornerHelper(negDirA, negDirB, dirC)
    addCornerHelper(negDirA, dirB, dirC)
    addCornerHelper(negDirA, dirB, negDirC)

    const matrix = Mat4.fromTranslation(Mat4(), Vec3.negate(Vec3(), origin))

    const cage = createCage(vertices, edges)
    const volume = Axes3D.volume(data.principalAxes.boxAxes)
    const radius = (Math.cbrt(volume) / 300) * props.scale
    state.currentGroup = VisualGroup.OrientedBox
    MeshBuilder.addCage(state, matrix, cage, radius, 2, 20)
}

function getOrientationMesh(data: OrientationData, props: OrientationProps, mesh?: Mesh) {
    const state = MeshBuilder.createState(256, 128, mesh)

    if (props.visuals.includes('principal-axes')) buildMomentsAxes(state, data, props)
    if (props.visuals.includes('oriented-box')) buildOrientedBox(state, data, props)

    return MeshBuilder.getMesh(state)
}

function getLabel(structure: Structure) {
    const loci = Structure.toStructureElementLoci(structure)
    const remappedLoci = StructureElement.Loci.remap(loci, structure.root)
    return structureElementStatsLabel(StructureElement.Stats.ofLoci(remappedLoci), { countsOnly: true })
}

export async function getStructureOrientationRepresentation(ctx: RuntimeContext, structure: Structure, params: OrientationProps, prev?: ShapeRepresentation<OrientationData, Mesh, Mesh.Params>) {
    const repr = prev || ShapeRepresentation(getOrientationShape, Mesh.Utils);
    const label = getLabel(structure)
    const principalAxes = PrincipalAxes.ofPositions(getPositions(structure))
    const data: OrientationData = { label, principalAxes }
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
    return Shape.create('Principal Axes', data, geo, () => props.color, () => 1, getLabel)
}