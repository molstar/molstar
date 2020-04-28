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
import { arcLength, halfPI, radToDeg } from '../../../mol-math/misc';
import { Circle } from '../../../mol-geo/primitive/circle';
import { transformPrimitive } from '../../../mol-geo/primitive/primitive';
import { MarkerActions, MarkerAction } from '../../../mol-util/marker-action';
import { dihedralLabel } from '../../../mol-theme/label';
import { MeasurementRepresentationCommonTextParams } from './common';
import { Sphere3D } from '../../../mol-math/geometry';

export interface DihedralData {
    quads: Loci.Bundle<4>[]
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

const ExtendersParams = {
    ...LinesParams
};
type ExtendersParams = typeof ExtendersParams

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
    borderWidth: PD.Numeric(0.2, { min: 0, max: 0.5, step: 0.01 }),
    ...MeasurementRepresentationCommonTextParams
};
type TextParams = typeof TextParams

const DihedralVisuals = {
    'vectors': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, VectorsParams>) => ShapeRepresentation(getVectorsShape, Lines.Utils, { modifyState: s => ({ ...s, pickable: false }) }),
    'extenders': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, ExtendersParams>) => ShapeRepresentation(getExtendersShape, Lines.Utils, { modifyState: s => ({ ...s, pickable: false }) }),
    'connector': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, ExtendersParams>) => ShapeRepresentation(getConnectorShape, Lines.Utils, { modifyState: s => ({ ...s, pickable: false }) }),
    'arc': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, ArcParams>) => ShapeRepresentation(getArcShape, Lines.Utils, { modifyState: s => ({ ...s, pickable: false }) }),
    'sector': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, SectorParams>) => ShapeRepresentation(getSectorShape, Mesh.Utils, { modifyProps: p => ({ ...p, alpha: p.sectorOpacity }), modifyState: s => ({ ...s, markerActions: MarkerActions.Highlighting }) }),
    'text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, TextParams>) => ShapeRepresentation(getTextShape, Text.Utils, { modifyState: s => ({ ...s, markerActions: MarkerAction.None }) }),
};

export const DihedralParams = {
    ...VectorsParams,
    ...ExtendersParams,
    ...ArcParams,
    ...SectorParams,
    ...TextParams,
    visuals: PD.MultiSelect(['extenders', 'sector', 'text'], PD.objectToOptions(DihedralVisuals)),
};
export type DihedralParams = typeof DihedralParams
export type DihedralProps = PD.Values<DihedralParams>

//

function getDihedralState() {
    return {
        sphereA: Sphere3D(),
        sphereB: Sphere3D(),
        sphereC: Sphere3D(),
        sphereD: Sphere3D(),

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
    };
}
type DihedralState = ReturnType<typeof getDihedralState>

const tmpVec = Vec3();
const tmpMat = Mat4();

// TODO improper dihedrals are not handled correctly
function setDihedralState(quad: Loci.Bundle<4>, state: DihedralState, arcScale: number) {
    const { sphereA, sphereB, sphereC, sphereD, dirBA, dirCD, projA, projD } = state;
    const { arcPointA, arcPointD, arcDirA, arcDirD, arcCenter, arcNormal } = state;

    const [lociA, lociB, lociC, lociD] = quad.loci;
    Loci.getBoundingSphere(lociA, sphereA);
    Loci.getBoundingSphere(lociB, sphereB);
    Loci.getBoundingSphere(lociC, sphereC);
    Loci.getBoundingSphere(lociD, sphereD);

    Vec3.add(arcCenter, sphereB.center, sphereC.center);
    Vec3.scale(arcCenter, arcCenter, 0.5);

    Vec3.sub(dirBA, sphereA.center, sphereB.center);
    Vec3.sub(dirCD, sphereD.center, sphereC.center);
    Vec3.add(arcPointA, arcCenter, dirBA);
    Vec3.add(arcPointD, arcCenter, dirCD);

    Vec3.sub(arcNormal, sphereC.center, sphereB.center);
    Vec3.orthogonalize(arcDirA, arcNormal, dirBA);
    Vec3.orthogonalize(arcDirD, arcNormal, dirCD);

    Vec3.projectPointOnVector(projA, arcPointA, arcDirA, arcCenter);
    Vec3.projectPointOnVector(projD, arcPointD, arcDirD, arcCenter);
    const len = Math.min(Vec3.distance(projA, arcCenter), Vec3.distance(projD, arcCenter));
    const radius = len * arcScale;

    Vec3.setMagnitude(arcDirA, arcDirA, radius);
    Vec3.setMagnitude(arcDirD, arcDirD, radius);
    Vec3.add(arcPointA, arcCenter, arcDirA);
    Vec3.add(arcPointD, arcCenter, arcDirD);
    state.radius = radius;

    state.angle = Vec3.dihedralAngle(sphereA.center, sphereB.center, sphereC.center, sphereD.center);

    Vec3.matchDirection(tmpVec, arcNormal, Vec3.sub(tmpVec, arcPointA, sphereA.center));
    const angleA = Vec3.angle(dirBA, tmpVec);
    const lenA = radius / Math.cos(angleA > halfPI ? angleA - halfPI : angleA);
    Vec3.add(projA, sphereB.center, Vec3.setMagnitude(tmpVec, dirBA, lenA));

    Vec3.matchDirection(tmpVec, arcNormal, Vec3.sub(tmpVec, arcPointD, sphereD.center));
    const angleD = Vec3.angle(dirCD, tmpVec);
    const lenD = radius / Math.cos(angleD > halfPI ? angleD - halfPI : angleD);
    Vec3.add(projD, sphereC.center, Vec3.setMagnitude(tmpVec, dirCD, lenD));

    return state;
}

function getCircle(state: DihedralState, segmentLength?: number) {
    const { radius, angle } = state;
    const segments = segmentLength ? arcLength(angle, radius) / segmentLength : 32;

    Mat4.targetTo(tmpMat, state.arcCenter, angle < 0 ? state.arcPointD : state.arcPointA, state.arcNormal);
    Mat4.setTranslation(tmpMat, state.arcCenter);
    Mat4.mul(tmpMat, tmpMat, Mat4.rotY180);

    const circle = Circle({ radius, thetaLength: Math.abs(angle), segments });
    return transformPrimitive(circle, tmpMat);
}

const tmpState = getDihedralState();

function getDihedralName(data: DihedralData) {
    return data.quads.length === 1 ? `Dihedral ${dihedralLabel(data.quads[0], { measureOnly: true })}` : `${data.quads.length} Dihedrals`;
}

//

function buildVectorsLines(data: DihedralData, props: DihedralProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines);
    for (let i = 0, il = data.quads.length; i < il; ++i) {
        setDihedralState(data.quads[i], tmpState, props.arcScale);
        builder.addFixedLengthDashes(tmpState.arcCenter, tmpState.arcPointA, props.dashLength, i);
        builder.addFixedLengthDashes(tmpState.arcCenter, tmpState.arcPointD, props.dashLength, i);
    }
    return builder.getLines();
}

function getVectorsShape(ctx: RuntimeContext, data: DihedralData, props: DihedralProps, shape?: Shape<Lines>) {
    const lines = buildVectorsLines(data, props, shape && shape.geometry);
    const name = getDihedralName(data);
    return Shape.create(name, data, lines, () => props.color, () => props.linesSize, () => '');
}

//

function buildConnectorLine(data: DihedralData, props: DihedralProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines);
    for (let i = 0, il = data.quads.length; i < il; ++i) {
        setDihedralState(data.quads[i], tmpState, props.arcScale);
        builder.addFixedLengthDashes(tmpState.sphereB.center, tmpState.sphereC.center, props.dashLength, i);
    }
    return builder.getLines();
}

function getConnectorShape(ctx: RuntimeContext, data: DihedralData, props: DihedralProps, shape?: Shape<Lines>) {
    const lines = buildConnectorLine(data, props, shape && shape.geometry);
    const name = getDihedralName(data);
    return Shape.create(name, data, lines, () => props.color, () => props.linesSize, () => '');
}

//

function buildExtendersLines(data: DihedralData, props: DihedralProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines);
    for (let i = 0, il = data.quads.length; i < il; ++i) {
        setDihedralState(data.quads[i], tmpState, props.arcScale);
        builder.addFixedLengthDashes(tmpState.arcPointA, tmpState.projA, props.dashLength, i);
        builder.addFixedLengthDashes(tmpState.arcPointD, tmpState.projD, props.dashLength, i);
    }
    return builder.getLines();
}

function getExtendersShape(ctx: RuntimeContext, data: DihedralData, props: DihedralProps, shape?: Shape<Lines>) {
    const lines = buildExtendersLines(data, props, shape && shape.geometry);
    const name = getDihedralName(data);
    return Shape.create(name, data, lines, () => props.color, () => props.linesSize, () => '');
}

//

function buildArcLines(data: DihedralData, props: DihedralProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines);
    for (let i = 0, il = data.quads.length; i < il; ++i) {
        setDihedralState(data.quads[i], tmpState, props.arcScale);
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

function getArcShape(ctx: RuntimeContext, data: DihedralData, props: DihedralProps, shape?: Shape<Lines>) {
    const lines = buildArcLines(data, props, shape && shape.geometry);
    const name = getDihedralName(data);
    return Shape.create(name, data, lines, () => props.color, () => props.linesSize, () => '');
}

//

function buildSectorMesh(data: DihedralData, props: DihedralProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(128, 64, mesh);
    for (let i = 0, il = data.quads.length; i < il; ++i) {
        setDihedralState(data.quads[i], tmpState, props.arcScale);
        const circle = getCircle(tmpState);
        state.currentGroup = i;
        MeshBuilder.addPrimitive(state, Mat4.id, circle);
        MeshBuilder.addPrimitiveFlipped(state, Mat4.id, circle);
    }
    return MeshBuilder.getMesh(state);
}

function getSectorShape(ctx: RuntimeContext, data: DihedralData, props: DihedralProps, shape?: Shape<Mesh>) {
    const mesh = buildSectorMesh(data, props, shape && shape.geometry);
    const name = getDihedralName(data);
    const getLabel = (groupId: number ) => dihedralLabel(data.quads[groupId]);
    return Shape.create(name, data, mesh, () => props.color, () => 1, getLabel);
}

//

function buildText(data: DihedralData, props: DihedralProps, text?: Text): Text {
    const builder = TextBuilder.create(props, 128, 64, text);
    for (let i = 0, il = data.quads.length; i < il; ++i) {
        setDihedralState(data.quads[i], tmpState, props.arcScale);

        Vec3.add(tmpVec, tmpState.arcDirA, tmpState.arcDirD);
        Vec3.setMagnitude(tmpVec, tmpVec, tmpState.radius);
        Vec3.add(tmpVec, tmpState.arcCenter, tmpVec);

        const angle = Math.abs(radToDeg(tmpState.angle)).toFixed(2);
        const label =  props.customText || `${angle}\u00B0`;
        const radius = Math.max(2, tmpState.sphereA.radius, tmpState.sphereB.radius, tmpState.sphereC.radius, tmpState.sphereD.radius);
        const scale = radius / 2;
        builder.add(label, tmpVec[0], tmpVec[1], tmpVec[2], 0.1, scale, i);
    }
    return builder.getText();
}

function getTextShape(ctx: RuntimeContext, data: DihedralData, props: DihedralProps, shape?: Shape<Text>) {
    const text = buildText(data, props, shape && shape.geometry);
    const name = getDihedralName(data);
    const getLabel = (groupId: number ) => dihedralLabel(data.quads[groupId]);
    return Shape.create(name, data, text, () => props.textColor, () => props.textSize, getLabel);
}

//

export type DihedralRepresentation = Representation<DihedralData, DihedralParams>
export function DihedralRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<DihedralData, DihedralParams>): DihedralRepresentation {
    return Representation.createMulti('Dihedral', ctx, getParams, Representation.StateBuilder, DihedralVisuals as unknown as Representation.Def<DihedralData, DihedralParams>);
}