/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../../../mol-model/loci';
import { RuntimeContext } from '../../../mol-task';
import { stringToWords } from '../../../mol-util/string';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { Text } from '../../../mol-geo/geometry/text/text';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ColorNames } from '../../../mol-util/color/names';
import { ShapeRepresentation } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../representation';
import { Shape } from '../../../mol-model/shape';
import { LinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { TextBuilder } from '../../../mol-geo/geometry/text/text-builder';
import { Vec3, Mat4 } from '../../../mol-math/linear-algebra';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { arcLength, halfPI, radToDeg } from '../../../mol-math/misc';
import { Circle } from '../../../mol-geo/primitive/circle';
import { transformPrimitive } from '../../../mol-geo/primitive/primitive';

export interface DihedralData {
    quads: Loci.Quad[]
}

const SharedParams = {
    color: PD.Color(ColorNames.lightgreen),
    arcScale: PD.Numeric(0.7, { min: 0.01, max: 1, step: 0.01 })
}

const LinesParams = {
    ...Lines.Params,
    ...SharedParams,
    lineSizeAttenuation: PD.Boolean(true),
    linesSize: PD.Numeric(0.04, { min: 0.01, max: 5, step: 0.01 }),
    dashLength: PD.Numeric(0.04, { min: 0.01, max: 0.2, step: 0.01 }),
}

const VectorsParams = {
    ...LinesParams
}
type VectorsParams = typeof VectorsParams

const ExtendersParams = {
    ...LinesParams
}
type ExtendersParams = typeof ExtendersParams

const ArcParams = {
    ...LinesParams
}
type ArcParams = typeof ArcParams

const SectorParams = {
    ...Mesh.Params,
    ...SharedParams,
    ignoreLight: PD.Boolean(true),
    sectorOpacity: PD.Numeric(0.75, { min: 0, max: 1, step: 0.01 }),
}
type SectorParams = typeof SectorParams

const TextParams = {
    ...Text.Params,
    borderWidth: PD.Numeric(0.2, { min: 0, max: 0.5, step: 0.01 }),
    textColor: PD.Color(ColorNames.black),
    textSize: PD.Numeric(0.4, { min: 0.1, max: 5, step: 0.1 }),
}
type TextParams = typeof TextParams

const DihedralVisuals = {
    'vectors': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, VectorsParams>) => ShapeRepresentation(getVectorsShape, Lines.Utils, { modifyState: s => ({ ...s, pickable: false }) }),
    'extenders': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, ExtendersParams>) => ShapeRepresentation(getExtendersShape, Lines.Utils, { modifyState: s => ({ ...s, pickable: false }) }),
    'arc': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, ArcParams>) => ShapeRepresentation(getArcShape, Lines.Utils, { modifyState: s => ({ ...s, pickable: false }) }),
    'sector': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, SectorParams>) => ShapeRepresentation(getSectorShape, Mesh.Utils, { modifyProps: p => ({ ...p, alpha: p.sectorOpacity }) }),
    'text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, TextParams>) => ShapeRepresentation(getTextShape, Text.Utils),
}
type DihedralVisualName = keyof typeof DihedralVisuals
const DihedralVisualOptions = Object.keys(DihedralVisuals).map(name => [name, stringToWords(name)] as [DihedralVisualName, string])

export const DihedralParams = {
    ...VectorsParams,
    ...ExtendersParams,
    ...ArcParams,
    ...SectorParams,
    ...TextParams,
    visuals: PD.MultiSelect<DihedralVisualName>(['extenders', 'sector', 'text'], DihedralVisualOptions),
}
export type DihedralParams = typeof DihedralParams
export type DihedralProps = PD.Values<DihedralParams>

//

function getDihedralState() {
    return {
        pointA: Vec3(),
        pointB: Vec3(),
        pointC: Vec3(),
        pointD: Vec3(),

        dirBA: Vec3(),
        dirCD: Vec3(),

        projA: Vec3(),
        projD: Vec3(),

        arcPointA: Vec3(),
        arcPointD: Vec3(),
        arcDirA: Vec3(),
        arcDirD: Vec3(),
        arcCenter: Vec3(),
        arcNormal: Vec3(),

        radius: 0,
        angle: 0,
    }
}
type DihedralState = ReturnType<typeof getDihedralState>

const tmpVec = Vec3()
const tmpMat = Mat4()

// TODO improper dihedrals are not handled correctly
function setDihedralState(quad: Loci.Quad, state: DihedralState, arcScale: number) {
    const { pointA, pointB, pointC, pointD, dirBA, dirCD, projA, projD } = state
    const { arcPointA, arcPointD, arcDirA, arcDirD, arcCenter, arcNormal } = state

    const { lociA, lociB, lociC, lociD } = quad
    Loci.getCenter(lociA, pointA)
    Loci.getCenter(lociB, pointB)
    Loci.getCenter(lociC, pointC)
    Loci.getCenter(lociD, pointD)

    Vec3.add(arcCenter, pointB, pointC)
    Vec3.scale(arcCenter, arcCenter, 0.5)

    Vec3.sub(dirBA, pointA, pointB)
    Vec3.sub(dirCD, pointD, pointC)
    Vec3.add(arcPointA, arcCenter, dirBA)
    Vec3.add(arcPointD, arcCenter, dirCD)

    Vec3.sub(arcNormal, pointC, pointB)
    Vec3.orthogonalize(arcDirA, arcNormal, dirBA)
    Vec3.orthogonalize(arcDirD, arcNormal, dirCD)

    Vec3.projectPointOnVector(projA, arcPointA, arcDirA, arcCenter)
    Vec3.projectPointOnVector(projD, arcPointD, arcDirD, arcCenter)
    const len = Math.min(Vec3.distance(projA, arcCenter), Vec3.distance(projD, arcCenter))
    const radius = len * arcScale

    Vec3.setMagnitude(arcDirA, arcDirA, radius)
    Vec3.setMagnitude(arcDirD, arcDirD, radius)
    Vec3.add(arcPointA, arcCenter, arcDirA)
    Vec3.add(arcPointD, arcCenter, arcDirD)
    state.radius = radius
    state.angle = Vec3.angle(arcDirA, arcDirD)

    Vec3.matchDirection(tmpVec, arcNormal, Vec3.sub(tmpVec, arcPointA, pointA))
    const angleA = Vec3.angle(dirBA, tmpVec)
    const lenA = radius / Math.cos(angleA > halfPI ? angleA - halfPI : angleA)
    Vec3.add(projA, pointB, Vec3.setMagnitude(tmpVec, dirBA, lenA))

    Vec3.matchDirection(tmpVec, arcNormal, Vec3.sub(tmpVec, arcPointD, pointD))
    const angleD = Vec3.angle(dirCD, tmpVec)
    const lenD = radius / Math.cos(angleD > halfPI ? angleD - halfPI : angleD)
    Vec3.add(projD, pointC, Vec3.setMagnitude(tmpVec, dirCD, lenD))

    return state
}

function getCircle(state: DihedralState, segmentLength?: number) {
    const { radius, angle } = state
    const segments = segmentLength ? arcLength(angle, radius) / segmentLength : 32

    Mat4.targetTo(tmpMat, state.arcCenter, angle > halfPI ? state.arcPointA : state.arcPointD, state.arcNormal)
    Mat4.setTranslation(tmpMat, state.arcCenter)
    Mat4.mul(tmpMat, tmpMat, Mat4.rotY180)

    const circle = Circle({ radius, thetaLength: angle, segments })
    return transformPrimitive(circle, tmpMat)
}

const tmpState = getDihedralState()
function dihedralLabel(quad: Loci.Quad, arcScale: number) {
    setDihedralState(quad, tmpState, arcScale)
    const angle = radToDeg(tmpState.angle).toFixed(2)
    return `Dihedral ${angle}\u00B0`
}

//

function buildVectorsLines(data: DihedralData, props: DihedralProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines)
    for (let i = 0, il = data.quads.length; i < il; ++i) {
        setDihedralState(data.quads[i], tmpState, props.arcScale)
        builder.addFixedLengthDashes(tmpState.arcCenter, tmpState.arcPointA, props.dashLength, i)
        builder.addFixedLengthDashes(tmpState.arcCenter, tmpState.arcPointD, props.dashLength, i)
    }
    return builder.getLines()
}

function getVectorsShape(ctx: RuntimeContext, data: DihedralData, props: DihedralProps, shape?: Shape<Lines>) {
    const lines = buildVectorsLines(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return dihedralLabel(data.quads[groupId], props.arcScale)
    }
    return Shape.create('Dihedral Vectors', data, lines, () => props.color, () => props.linesSize, getLabel)
}

//

function buildExtendersLines(data: DihedralData, props: DihedralProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines)
    for (let i = 0, il = data.quads.length; i < il; ++i) {
        setDihedralState(data.quads[i], tmpState, props.arcScale)
        builder.addFixedLengthDashes(tmpState.arcPointA, tmpState.projA, props.dashLength, i)
        builder.addFixedLengthDashes(tmpState.arcPointD, tmpState.projD, props.dashLength, i)
    }
    return builder.getLines()
}

function getExtendersShape(ctx: RuntimeContext, data: DihedralData, props: DihedralProps, shape?: Shape<Lines>) {
    const lines = buildExtendersLines(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return dihedralLabel(data.quads[groupId], props.arcScale)
    }
    return Shape.create('Dihedral Extenders', data, lines, () => props.color, () => props.linesSize, getLabel)
}

//

function buildArcLines(data: DihedralData, props: DihedralProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines)
    for (let i = 0, il = data.quads.length; i < il; ++i) {
        setDihedralState(data.quads[i], tmpState, props.arcScale)
        const circle = getCircle(tmpState, props.dashLength)
        const { indices, vertices } = circle
        for (let j = 0, jl = indices.length; j < jl; j += 3) {
            if (j % 2 === 1) continue // draw every other segment to get dashes
            const start = indices[j] * 3
            const end = indices[j + 1] * 3
            const startX = vertices[start]
            const startY = vertices[start + 1]
            const startZ = vertices[start + 2]
            const endX = vertices[end]
            const endY = vertices[end + 1]
            const endZ = vertices[end + 2]
            builder.add(startX, startY, startZ, endX, endY, endZ, i)
        }
    }
    return builder.getLines()
}

function getArcShape(ctx: RuntimeContext, data: DihedralData, props: DihedralProps, shape?: Shape<Lines>) {
    const lines = buildArcLines(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return dihedralLabel(data.quads[groupId], props.arcScale)
    }
    return Shape.create('Dihedral Arc', data, lines, () => props.color, () => props.linesSize, getLabel)
}

//

function buildSectorMesh(data: DihedralData, props: DihedralProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(128, 64, mesh)
    for (let i = 0, il = data.quads.length; i < il; ++i) {
        setDihedralState(data.quads[i], tmpState, props.arcScale)
        const circle = getCircle(tmpState)
        MeshBuilder.addPrimitive(state, Mat4.id, circle)
        MeshBuilder.addPrimitiveFlipped(state, Mat4.id, circle)
    }
    return MeshBuilder.getMesh(state)
}

function getSectorShape(ctx: RuntimeContext, data: DihedralData, props: DihedralProps, shape?: Shape<Mesh>) {
    const mesh = buildSectorMesh(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return dihedralLabel(data.quads[groupId], props.arcScale)
    }
    return Shape.create('Dihedral Sector', data, mesh, () => props.color, () => 1, getLabel)
}

//

function buildText(data: DihedralData, props: DihedralProps, text?: Text): Text {
    const builder = TextBuilder.create(props, 128, 64, text)
    for (let i = 0, il = data.quads.length; i < il; ++i) {
        setDihedralState(data.quads[i], tmpState, props.arcScale)

        Vec3.add(tmpVec, tmpState.arcDirA, tmpState.arcDirD)
        Vec3.setMagnitude(tmpVec, tmpVec, tmpState.radius)
        Vec3.add(tmpVec, tmpState.arcCenter, tmpVec)

        const angle = radToDeg(tmpState.angle).toFixed(2)
        const label = `${angle}\u00B0`
        builder.add(label, tmpVec[0], tmpVec[1], tmpVec[2], 0.1, 1, i)
    }
    return builder.getText()
}

function getTextShape(ctx: RuntimeContext, data: DihedralData, props: DihedralProps, shape?: Shape<Text>) {
    const text = buildText(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return dihedralLabel(data.quads[groupId], props.arcScale)
    }
    return Shape.create('Dihedral Text', data, text, () => props.textColor, () => props.textSize, getLabel)
}

//

export type DihedralRepresentation = Representation<DihedralData, DihedralParams>
export function DihedralRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, DihedralParams>): DihedralRepresentation {
    return Representation.createMulti('Dihedral', ctx, getParams, Representation.StateBuilder, DihedralVisuals as unknown as Representation.Def<DihedralData, DihedralParams>)
}