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
import { radToDeg } from '../../../mol-math/misc';
import { Circle } from '../../../mol-geo/primitive/circle';

export interface AngleData {
    triplets: { lociA: Loci, lociB: Loci, lociC: Loci }[]
}

const SharedParams = {
    color: PD.Color(ColorNames.darkgreen),
}

const LinesParams = {
    ...Lines.Params,
    ...SharedParams,
    lineSizeAttenuation: PD.Boolean(true),
    linesSize: PD.Numeric(0.05, { min: 0.01, max: 5, step: 0.01 }),
}

const VectorsParams = {
    ...LinesParams
}
type VectorsParams = typeof VectorsParams

const ArcParams = {
    ...LinesParams
}
type ArcParams = typeof ArcParams

const SectorParams = {
    ...Mesh.Params,
    ...SharedParams,
}
type SectorParams = typeof SectorParams

const TextParams = {
    ...Text.Params,
    borderWidth: PD.Numeric(0.25, { min: 0, max: 0.5, step: 0.01 }),
    textColor: PD.Color(ColorNames.black),
    textSize: PD.Numeric(0.8, { min: 0.1, max: 5, step: 0.1 }),
}
type TextParams = typeof TextParams

const AngleVisuals = {
    'vectors': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<AngleData, VectorsParams>) => ShapeRepresentation(getVectorsShape, Lines.Utils),
    'arc': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<AngleData, ArcParams>) => ShapeRepresentation(getArcShape, Lines.Utils),
    'sector': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<AngleData, SectorParams>) => ShapeRepresentation(getSectorShape, Mesh.Utils),
    'text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<AngleData, TextParams>) => ShapeRepresentation(getTextShape, Text.Utils),
}
type AngleVisualName = keyof typeof AngleVisuals
const AngleVisualOptions = Object.keys(AngleVisuals).map(name => [name, stringToWords(name)] as [AngleVisualName, string])

export const AngleParams = {
    ...VectorsParams,
    ...ArcParams,
    ...SectorParams,
    ...TextParams,
    visuals: PD.MultiSelect<AngleVisualName>(['vectors', 'arc', 'sector', 'text'], AngleVisualOptions),
}
export type AngleParams = typeof AngleParams
export type AngleProps = PD.Values<AngleParams>

//

const tmpVecA = Vec3()
const tmpVecB = Vec3()
const tmpVecC = Vec3()
// const tmpVecD = Vec3()

const tmpDirA = Vec3()
const tmpDirB = Vec3()
const tmpCenter = Vec3()

//

function buildVectorsLines(data: AngleData, props: AngleProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines)
    for (let i = 0, il = data.triplets.length; i < il; ++i) {
        const { lociA, lociB, lociC } = data.triplets[i]
        Loci.getCenter(lociA, tmpVecA)
        Loci.getCenter(lociB, tmpVecB)
        Loci.getCenter(lociC, tmpVecC)
        builder.addFixedLengthDashes(tmpVecA, tmpVecB, 0.1, i)
        builder.addFixedLengthDashes(tmpVecB, tmpVecC, 0.1, i)
    }
    return builder.getLines()
}

function getVectorsShape(ctx: RuntimeContext, data: AngleData, props: AngleProps, shape?: Shape<Lines>) {
    const lines = buildVectorsLines(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return 'Angle Vectors'
    }
    return Shape.create('Angle Vectors', data, lines, () => props.color, () => props.linesSize, getLabel)
}

//

function buildArcLines(data: AngleData, props: AngleProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines)

    return builder.getLines()
}

function getArcShape(ctx: RuntimeContext, data: AngleData, props: AngleProps, shape?: Shape<Lines>) {
    const lines = buildArcLines(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return 'Angle Arc'
    }
    return Shape.create('Angle Arc', data, lines, () => props.color, () => props.linesSize, getLabel)
}

//

function buildSectorMesh(data: AngleData, props: AngleProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(128, 64, mesh)
    const m = Mat4()
    const tmpVec = Vec3()

    for (let i = 0, il = data.triplets.length; i < il; ++i) {
        const { lociA, lociB, lociC } = data.triplets[i]
        console.log(data.triplets[i])
        Loci.getCenter(lociA, tmpVecA)
        Loci.getCenter(lociB, tmpVecB)
        Loci.getCenter(lociC, tmpVecC)
        Vec3.sub(tmpDirA, tmpVecA, tmpVecB)
        Vec3.sub(tmpDirB, tmpVecC, tmpVecB)

        const lenA = Vec3.magnitude(tmpDirA)
        const lenB = Vec3.magnitude(tmpDirB)
        let dir: Vec3, len: number
        if (lenA <= lenB) {
            dir = tmpDirA
            len = lenA
        } else {
            dir = tmpDirB
            len = lenB
        }

        const center = tmpVecB
        const dirMajor = Vec3.cross(Vec3(), tmpDirA, tmpDirB)
        const dirMinor = dir

        Vec3.add(tmpVec, center, dirMajor)
        Mat4.targetTo(m, center, tmpVec, dirMinor)
        Mat4.setTranslation(m, center)

        const angle = Vec3.angle(tmpDirA, tmpDirB)
        const circle = Circle({ radius: len, thetaLength: angle })
        console.log(circle)

        MeshBuilder.addPrimitive(state, m, circle)
        MeshBuilder.addPrimitiveFlipped(state, m, circle)
    }

    return MeshBuilder.getMesh(state)
}

function getSectorShape(ctx: RuntimeContext, data: AngleData, props: AngleProps, shape?: Shape<Mesh>) {
    const mesh = buildSectorMesh(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return 'Angle Sector'
    }
    return Shape.create('Angle Sector', data, mesh, () => props.color, () => 1, getLabel)
}

//

function buildText(data: AngleData, props: AngleProps, text?: Text): Text {
    const builder = TextBuilder.create(props, 128, 64, text)
    for (let i = 0, il = data.triplets.length; i < il; ++i) {
        const {  lociA, lociB, lociC } = data.triplets[i]
        Loci.getCenter(lociA, tmpVecA)
        Loci.getCenter(lociB, tmpVecB)
        Loci.getCenter(lociC, tmpVecC)
        Vec3.sub(tmpDirA, tmpVecB, tmpVecA)
        Vec3.sub(tmpDirB, tmpVecC, tmpVecB)
        Vec3.add(tmpCenter, tmpVecA, tmpVecC)
        Vec3.scale(tmpCenter, tmpCenter, 0.5)
        const angle = radToDeg(Vec3.angle(tmpDirA, tmpDirB)).toPrecision(2)
        const label = `${angle}\u00B0`
        builder.add(label, tmpCenter[0], tmpCenter[1], tmpCenter[2], 0.1, i)
    }
    return builder.getText()
}

function getTextShape(ctx: RuntimeContext, data: AngleData, props: AngleProps, shape?: Shape<Text>) {
    const text = buildText(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return 'Angle Text'
    }
    return Shape.create('Angle Text', data, text, () => props.textColor, () => props.textSize, getLabel)
}

//

export type AngleRepresentation = Representation<AngleData, AngleParams>
export function AngleRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<AngleData, AngleParams>): AngleRepresentation {
    return Representation.createMulti('Angle', ctx, getParams, Representation.StateBuilder, AngleVisuals as unknown as Representation.Def<AngleData, AngleParams>)
}