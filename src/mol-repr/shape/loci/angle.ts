/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../../../mol-model/loci';
import { RuntimeContext } from '../../../mol-task';
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
import { radToDeg, arcLength } from '../../../mol-math/misc';
import { Circle } from '../../../mol-geo/primitive/circle';
import { transformPrimitive } from '../../../mol-geo/primitive/primitive';
import { MarkerActions, MarkerAction } from '../../../mol-util/marker-action';
import { angleLabel } from '../../../mol-theme/label';
import { Sphere3D } from '../../../mol-math/geometry';
import { MeasurementRepresentationCommonTextParams } from './common';

export interface AngleData {
    triples: Loci.Bundle<3>[]
}

const SharedParams = {
    color: PD.Color(ColorNames.lightgreen),
    arcScale: PD.Numeric(0.7, { min: 0.01, max: 1, step: 0.01 })
};

const LinesParams = {
    ...Lines.Params,
    ...SharedParams,
    lineSizeAttenuation: PD.Boolean(true),
    linesSize: PD.Numeric(0.04, { min: 0.01, max: 5, step: 0.01 }),
    dashLength: PD.Numeric(0.04, { min: 0.01, max: 0.2, step: 0.01 }),
};

const VectorsParams = {
    ...LinesParams
};
type VectorsParams = typeof VectorsParams

const ArcParams = {
    ...LinesParams
};
type ArcParams = typeof ArcParams

const SectorParams = {
    ...Mesh.Params,
    ...SharedParams,
    ignoreLight: PD.Boolean(true),
    sectorOpacity: PD.Numeric(0.75, { min: 0, max: 1, step: 0.01 }),
};
type SectorParams = typeof SectorParams

const TextParams = {
    ...Text.Params,
    ...MeasurementRepresentationCommonTextParams,
    borderWidth: PD.Numeric(0.2, { min: 0, max: 0.5, step: 0.01 })
};
type TextParams = typeof TextParams

const AngleVisuals = {
    'vectors': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<AngleData, VectorsParams>) => ShapeRepresentation(getVectorsShape, Lines.Utils, { modifyState: s => ({ ...s, pickable: false }) }),
    'arc': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<AngleData, ArcParams>) => ShapeRepresentation(getArcShape, Lines.Utils, { modifyState: s => ({ ...s, pickable: false }) }),
    'sector': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<AngleData, SectorParams>) => ShapeRepresentation(getSectorShape, Mesh.Utils, { modifyProps: p => ({ ...p, alpha: p.sectorOpacity }), modifyState: s => ({ ...s, markerActions: MarkerActions.Highlighting }) }),
    'text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<AngleData, TextParams>) => ShapeRepresentation(getTextShape, Text.Utils, { modifyState: s => ({ ...s, markerActions: MarkerAction.None }) }),
};

export const AngleParams = {
    ...VectorsParams,
    ...ArcParams,
    ...SectorParams,
    ...TextParams,
    visuals: PD.MultiSelect(['vectors', 'sector', 'text'], PD.objectToOptions(AngleVisuals)),
};
export type AngleParams = typeof AngleParams
export type AngleProps = PD.Values<AngleParams>

//

function getAngleState() {
    return {
        sphereA: Sphere3D(),
        sphereB: Sphere3D(),
        sphereC: Sphere3D(),

        arcDirA: Vec3(),
        arcDirC: Vec3(),
        arcNormal: Vec3(),

        radius: 0,
        angle: 0,
    };
}
type AngleState = ReturnType<typeof getAngleState>

const tmpVec = Vec3();
const tmpMat = Mat4();

function setAngleState(triple: Loci.Bundle<3>, state: AngleState, arcScale: number) {
    const { sphereA, sphereB, sphereC } = state;
    const { arcDirA, arcDirC, arcNormal } = state;

    const [lociA, lociB, lociC] = triple.loci;
    Loci.getBoundingSphere(lociA, sphereA);
    Loci.getBoundingSphere(lociB, sphereB);
    Loci.getBoundingSphere(lociC, sphereC);

    Vec3.sub(arcDirA, sphereA.center, sphereB.center);
    Vec3.sub(arcDirC, sphereC.center, sphereB.center);
    Vec3.cross(arcNormal, arcDirA, arcDirC);

    const len = Math.min(Vec3.magnitude(arcDirA), Vec3.magnitude(arcDirC));
    const radius = len * arcScale;

    state.radius = radius;
    state.angle = Vec3.angle(arcDirA, arcDirC);

    return state;
}

function getCircle(state: AngleState, segmentLength?: number) {
    const { radius, angle } = state;
    const segments = segmentLength ? arcLength(angle, radius) / segmentLength : 32;

    Mat4.targetTo(tmpMat, state.sphereB.center, state.sphereA.center, state.arcNormal);
    Mat4.setTranslation(tmpMat, state.sphereB.center);
    Mat4.mul(tmpMat, tmpMat, Mat4.rotY180);

    const circle = Circle({ radius, thetaLength: angle, segments });
    return transformPrimitive(circle, tmpMat);
}

const tmpState = getAngleState();

function getAngleName(data: AngleData) {
    return data.triples.length === 1 ? `Angle ${angleLabel(data.triples[0], { measureOnly: true })}` : `${data.triples.length} Angles`;
}

//

function buildVectorsLines(data: AngleData, props: AngleProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines);
    for (let i = 0, il = data.triples.length; i < il; ++i) {
        setAngleState(data.triples[i], tmpState, props.arcScale);
        builder.addFixedLengthDashes(tmpState.sphereB.center, tmpState.sphereA.center, props.dashLength, i);
        builder.addFixedLengthDashes(tmpState.sphereB.center, tmpState.sphereC.center, props.dashLength, i);
    }
    return builder.getLines();
}

function getVectorsShape(ctx: RuntimeContext, data: AngleData, props: AngleProps, shape?: Shape<Lines>) {
    const lines = buildVectorsLines(data, props, shape && shape.geometry);
    const name = getAngleName(data);
    return Shape.create(name, data, lines, () => props.color, () => props.linesSize, () => '');
}

//

function buildArcLines(data: AngleData, props: AngleProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines);
    for (let i = 0, il = data.triples.length; i < il; ++i) {
        setAngleState(data.triples[i], tmpState, props.arcScale);
        const circle = getCircle(tmpState, props.dashLength);
        const { indices, vertices } = circle;
        for (let j = 0, jl = indices.length; j < jl; j += 3) {
            if (j % 2 === 1) continue; // draw every other segment to get dashes
            const start = indices[j] * 3;
            const end = indices[j + 1] * 3;
            const startX = vertices[start];
            const startY = vertices[start + 1];
            const startZ = vertices[start + 2];
            const endX = vertices[end];
            const endY = vertices[end + 1];
            const endZ = vertices[end + 2];
            builder.add(startX, startY, startZ, endX, endY, endZ, i);
        }
    }
    return builder.getLines();
}

function getArcShape(ctx: RuntimeContext, data: AngleData, props: AngleProps, shape?: Shape<Lines>) {
    const lines = buildArcLines(data, props, shape && shape.geometry);
    const name = getAngleName(data);
    return Shape.create(name, data, lines, () => props.color, () => props.linesSize, () => '');
}

//

function buildSectorMesh(data: AngleData, props: AngleProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(128, 64, mesh);
    for (let i = 0, il = data.triples.length; i < il; ++i) {
        setAngleState(data.triples[i], tmpState, props.arcScale);
        const circle = getCircle(tmpState);
        state.currentGroup = i;
        MeshBuilder.addPrimitive(state, Mat4.id, circle);
        MeshBuilder.addPrimitiveFlipped(state, Mat4.id, circle);
    }
    return MeshBuilder.getMesh(state);
}

function getSectorShape(ctx: RuntimeContext, data: AngleData, props: AngleProps, shape?: Shape<Mesh>) {
    const mesh = buildSectorMesh(data, props, shape && shape.geometry);
    const name = getAngleName(data);
    const getLabel = (groupId: number ) => angleLabel(data.triples[groupId]);
    return Shape.create(name, data, mesh, () => props.color, () => 1, getLabel);
}

//

function buildText(data: AngleData, props: AngleProps, text?: Text): Text {
    const builder = TextBuilder.create(props, 128, 64, text);
    for (let i = 0, il = data.triples.length; i < il; ++i) {
        setAngleState(data.triples[i], tmpState, props.arcScale);

        Vec3.add(tmpVec, tmpState.arcDirA, tmpState.arcDirC);
        Vec3.setMagnitude(tmpVec, tmpVec, tmpState.radius);
        Vec3.add(tmpVec, tmpState.sphereB.center, tmpVec);

        const angle = radToDeg(tmpState.angle).toFixed(2);
        const label = props.customText || `${angle}\u00B0`;
        const radius = Math.max(2, tmpState.sphereA.radius, tmpState.sphereB.radius, tmpState.sphereC.radius);
        const scale = radius / 2;
        builder.add(label, tmpVec[0], tmpVec[1], tmpVec[2], 0.1, scale, i);
    }
    return builder.getText();
}

function getTextShape(ctx: RuntimeContext, data: AngleData, props: AngleProps, shape?: Shape<Text>) {
    const text = buildText(data, props, shape && shape.geometry);
    const name = getAngleName(data);
    const getLabel = (groupId: number ) => angleLabel(data.triples[groupId]);
    return Shape.create(name, data, text, () => props.textColor, () => props.textSize, getLabel);
}

//

export type AngleRepresentation = Representation<AngleData, AngleParams>
export function AngleRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<AngleData, AngleParams>): AngleRepresentation {
    return Representation.createMulti('Angle', ctx, getParams, Representation.StateBuilder, AngleVisuals as unknown as Representation.Def<AngleData, AngleParams>);
}