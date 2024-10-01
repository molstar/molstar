/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Vec3, Mat4 } from '../../../../mol-math/linear-algebra';
import { MeshBuilder } from '../mesh-builder';
import { Primitive, transformPrimitive } from '../../../primitive/primitive';
import { Line, LineProps, DefaultLineProps } from '../../../primitive/line';
import { Prism } from '../../../primitive/prism';
import { polygon } from '../../../primitive/polygon';
import { hashFnv32a } from '../../../../mol-data/util';

const lineMap = new Map<number, Primitive>();
const up = Vec3.create(0, 1, 0);

const tmpLineDir = Vec3();
const tmpLineMatDir = Vec3();
const tmpLineCenter = Vec3();
const tmpLineMat = Mat4();
const tmpLineMatRot = Mat4();
const tmpLineScale = Vec3();
const tmpLineStart = Vec3();
const tmpUp = Vec3();

function setLineMat(m: Mat4, start: Vec3, dir: Vec3, length: number, matchDir: boolean) {
    Vec3.setMagnitude(tmpLineMatDir, dir, length / 2);
    Vec3.add(tmpLineCenter, start, tmpLineMatDir);
    // ensure the direction used to create the rotation is always pointing in the same
    // direction so the triangles of adjacent cylinder will line up
    if (matchDir) Vec3.matchDirection(tmpUp, up, tmpLineMatDir);
    else Vec3.copy(tmpUp, up);
    Vec3.set(tmpLineScale, 1, length, 1);
    Vec3.makeRotation(tmpLineMatRot, tmpUp, tmpLineMatDir);
    Mat4.scale(m, tmpLineMatRot, tmpLineScale);
    return Mat4.setTranslation(m, tmpLineCenter);
}

const tmpPropValues = new Int32Array(9);
function getLinePropsKey(props: LineProps) {
    // TODO: remove
    const radiusTop = 1;
    const radiusBottom = 1;
    const height = 1;
    tmpPropValues[0] = Math.round(1000 * (radiusTop));
    tmpPropValues[1] = Math.round(1000 * (radiusBottom));
    tmpPropValues[2] = Math.round(1000 * (height));
    tmpPropValues[3] = props.radialSegments ?? DefaultLineProps.radialSegments;
    tmpPropValues[4] = props.heightSegments ?? DefaultLineProps.heightSegments;
    tmpPropValues[5] = (props.topCap ?? DefaultLineProps.topCap) ? 1 : 0;
    tmpPropValues[6] = (props.bottomCap ?? DefaultLineProps.bottomCap) ? 1 : 0;
    tmpPropValues[7] = Math.round(1000 * (props.thetaStart ?? DefaultLineProps.thetaStart));
    tmpPropValues[8] = Math.round(1000 * (props.thetaLength ?? DefaultLineProps.thetaLength));
    return hashFnv32a(tmpPropValues);
}

function getLine(props: LineProps) {
    const radiusTop = 1;
    const key = getLinePropsKey(props);
    let line = lineMap.get(key);
    if (line === undefined) {
        if (props.radialSegments && props.radialSegments <= 4) {
            const sideCount = Math.max(3, props.radialSegments);
            const prism = Prism(polygon(sideCount, true, radiusTop), props);
            line = transformPrimitive(prism, Mat4.rotX90);
        } else {
            line = Line(props);
        }
        lineMap.set(key, line);
    }
    return line;
}

export type BasicLineProps = Omit<LineProps, 'height'>
export type StrokeDasharray = number | number[];

export function addLine(state: MeshBuilder.State, start: Vec3, end: Vec3, lengthScale: number, props: BasicLineProps, strokeDasharray?: StrokeDasharray) {
    const d = Vec3.distance(start, end) * lengthScale;
    Vec3.sub(tmpLineDir, end, start);
    setLineMat(tmpLineMat, start, tmpLineDir, d, true);
    MeshBuilder.addPrimitive(state, tmpLineMat, getLine(props));
}

export function addFixedCountDashedLine(state: MeshBuilder.State, start: Vec3, end: Vec3, lengthScale: number, segmentCount: number, stubCap: boolean, props: BasicLineProps) {
    const d = Vec3.distance(start, end) * lengthScale;
    const isOdd = segmentCount % 2 !== 0;
    const s = Math.floor((segmentCount + 1) / 2);
    let step = d / (segmentCount + 0.5);

    let line = getLine(props);
    Vec3.setMagnitude(tmpLineDir, Vec3.sub(tmpLineDir, end, start), step);
    Vec3.copy(tmpLineStart, start);
    for (let j = 0; j < s; ++j) {
        Vec3.add(tmpLineStart, tmpLineStart, tmpLineDir);
        if (isOdd && j === s - 1) {
            if (!stubCap && props.topCap) {
                props.topCap = false;
                line = getLine(props);
            }
            step /= 2;
        }
        setLineMat(tmpLineMat, tmpLineStart, tmpLineDir, step, false);
        MeshBuilder.addPrimitive(state, tmpLineMat, line);

        Vec3.add(tmpLineStart, tmpLineStart, tmpLineDir);
    }
}
