/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from '../../../../mol-math/linear-algebra';
import { MeshBuilder } from '../mesh-builder';
import { Primitive, transformPrimitive } from '../../../primitive/primitive';
import { Cylinder, CylinderProps } from '../../../primitive/cylinder';
import { Prism } from '../../../primitive/prism';
import { polygon } from '../../../primitive/polygon';

const cylinderMap = new Map<string, Primitive>();
const up = Vec3.create(0, 1, 0);

const tmpCylinderDir = Vec3.zero();
const tmpCylinderMatDir = Vec3.zero();
const tmpCylinderCenter = Vec3.zero();
const tmpCylinderMat = Mat4.zero();
const tmpCylinderStart = Vec3.zero();
const tmpUp = Vec3.zero();

function setCylinderMat(m: Mat4, start: Vec3, dir: Vec3, length: number) {
    Vec3.setMagnitude(tmpCylinderMatDir, dir, length / 2);
    Vec3.add(tmpCylinderCenter, start, tmpCylinderMatDir);
    // ensure the direction used to create the rotation is always pointing in the same
    // direction so the triangles of adjacent cylinder will line up
    Vec3.matchDirection(tmpUp, up, tmpCylinderMatDir);
    Vec3.makeRotation(m, tmpUp, tmpCylinderMatDir);
    return Mat4.setTranslation(m, tmpCylinderCenter);
}

function getCylinder(props: CylinderProps) {
    const key = JSON.stringify(props);
    let cylinder = cylinderMap.get(key);
    if (cylinder === undefined) {
        if (props.radialSegments && props.radialSegments <= 4) {
            const box = Prism(polygon(4, true, props.radiusTop), props);
            cylinder = transformPrimitive(box, Mat4.rotX90);
        } else {
            cylinder = Cylinder(props);
        }
        cylinderMap.set(key, cylinder);
    }
    return cylinder;
}

export function addCylinder(state: MeshBuilder.State, start: Vec3, end: Vec3, lengthScale: number, props: CylinderProps) {
    const d = Vec3.distance(start, end) * lengthScale;
    props.height = d;
    Vec3.sub(tmpCylinderDir, end, start);
    setCylinderMat(tmpCylinderMat, start, tmpCylinderDir, d);
    MeshBuilder.addPrimitive(state, tmpCylinderMat, getCylinder(props));
}

export function addDoubleCylinder(state: MeshBuilder.State, start: Vec3, end: Vec3, lengthScale: number, shift: Vec3, props: CylinderProps) {
    const d = Vec3.distance(start, end) * lengthScale;
    props.height = d;
    const cylinder = getCylinder(props);
    Vec3.sub(tmpCylinderDir, end, start);
    // positivly shifted cylinder
    Vec3.add(tmpCylinderStart, start, shift);
    setCylinderMat(tmpCylinderMat, tmpCylinderStart, tmpCylinderDir, d);
    MeshBuilder.addPrimitive(state, tmpCylinderMat, cylinder);
    // negativly shifted cylinder
    Vec3.sub(tmpCylinderStart, start, shift);
    setCylinderMat(tmpCylinderMat, tmpCylinderStart, tmpCylinderDir, d);
    MeshBuilder.addPrimitive(state, tmpCylinderMat, cylinder);
}

export function addFixedCountDashedCylinder(state: MeshBuilder.State, start: Vec3, end: Vec3, lengthScale: number, segmentCount: number, props: CylinderProps) {
    const s = Math.floor(segmentCount / 2);
    const step = 1 / segmentCount;

    // automatically adjust length so links/bonds that are rendered as two half cylinders
    // have evenly spaced dashed cylinders
    if (lengthScale < 1) {
        const bias = lengthScale / 2 / segmentCount;
        lengthScale += segmentCount % 2 === 1 ? bias : -bias;
    }

    const d = Vec3.distance(start, end) * lengthScale;
    props.height = d * step;
    const cylinder = getCylinder(props);
    Vec3.sub(tmpCylinderDir, end, start);

    for (let j = 0; j < s; ++j) {
        const f = step * (j * 2 + 1);
        Vec3.setMagnitude(tmpCylinderDir, tmpCylinderDir, d * f);
        Vec3.add(tmpCylinderStart, start, tmpCylinderDir);
        setCylinderMat(tmpCylinderMat, tmpCylinderStart, tmpCylinderDir, d * step);
        MeshBuilder.addPrimitive(state, tmpCylinderMat, cylinder);
    }
}