/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Quat, Mat4 } from '../../mol-math/linear-algebra';
import { NumberArray } from '../../mol-util/type-helpers';

interface Frame {
    t: Vec3,
    r: Vec3,
    s: Vec3,
}

const a0Tmp = Vec3();
const a1Tmp = Vec3();
const a2Tmp = Vec3();
const a3Tmp = Vec3();
function CubicInterpolate(out: Vec3, y0: Vec3, y1: Vec3, y2: Vec3, y3: Vec3, mu: number): Vec3 {
    const mu2 = mu * mu;
    Vec3.sub(a0Tmp, y3, y2);
    Vec3.sub(a0Tmp, a0Tmp, y0);
    Vec3.add(a0Tmp, a0Tmp, y1);

    Vec3.sub(a1Tmp, y0, y1);
    Vec3.sub(a1Tmp, a1Tmp, a0Tmp);

    Vec3.sub(a2Tmp, y2, y0);

    Vec3.copy(a3Tmp, y1);

    out[0] = a0Tmp[0] * mu * mu2 + a1Tmp[0] * mu2 + a2Tmp[0] * mu + a3Tmp[0];
    out[1] = a0Tmp[1] * mu * mu2 + a1Tmp[1] * mu2 + a2Tmp[1] * mu + a3Tmp[1];
    out[2] = a0Tmp[2] * mu * mu2 + a1Tmp[2] * mu2 + a2Tmp[2] * mu + a3Tmp[2];

    return out;
}

const cp0 = Vec3();
const cp1 = Vec3();
const cp2 = Vec3();
const cp3 = Vec3();
const currentPosition = Vec3();
function ResampleControlPoints(points: NumberArray, segmentLength: number) {
    const nP = points.length / 3;
    // insert a point at the end and at the begining
    // controlPoints.Insert(0, controlPoints[0] + (controlPoints[0] - controlPoints[1]) / 2.0f);
    // controlPoints.Add(controlPoints[nP - 1] + (controlPoints[nP - 1] - controlPoints[nP - 2]) / 2.0f);

    let resampledControlPoints: Vec3[] = [];
    // resampledControlPoints.Add(controlPoints[0]);
    // resampledControlPoints.Add(controlPoints[1]);

    let idx = 1;
    // const currentPosition = Vec3.create(points[idx * 3], points[idx * 3 + 1], points[idx * 3 + 2])
    Vec3.fromArray(currentPosition, points, idx * 3);

    let lerpValue = 0.0;

    // Normalize the distance between control points
    while (true) {
        if (idx + 2 >= nP) break;
        Vec3.fromArray(cp0, points, (idx - 1) * 3);
        Vec3.fromArray(cp1, points, idx * 3);
        Vec3.fromArray(cp2, points, (idx + 1) * 3);
        Vec3.fromArray(cp3, points, (idx + 2) * 3);
        // const cp0 = Vec3.create(points[(idx-1)*3], points[(idx-1)*3+1], points[(idx-1)*3+2]) // controlPoints[currentPointId - 1];
        // const cp1 = Vec3.create(points[idx*3], points[idx*3+1], points[idx*3+2]) // controlPoints[currentPointId];
        // const cp2 = Vec3.create(points[(idx+1)*3], points[(idx+1)*3+1], points[(idx+1)*3+2]) // controlPoints[currentPointId + 1];
        // const cp3 = Vec3.create(points[(idx+2)*3], points[(idx+2)*3+1], points[(idx+2)*3+2]); // controlPoints[currentPointId + 2];
        let found = false;
        for (; lerpValue <= 1; lerpValue += 0.01) {
            // lerp?slerp
            // let candidate:Vec3 = Vec3.lerp(Vec3.zero(), cp0, cp1, lerpValue);
            // const candidate:Vec3 = Vec3.bezier(Vec3.zero(), cp0, cp1, cp2, cp3, lerpValue);
            const candidate = CubicInterpolate(Vec3(), cp0, cp1, cp2, cp3, lerpValue);
            const d = Vec3.distance(currentPosition, candidate);
            if (d > segmentLength) {
                resampledControlPoints.push(candidate);
                Vec3.copy(currentPosition, candidate);
                found = true;
                break;
            }
        }
        if (!found) {
            lerpValue = 0;
            idx += 1;
        }
    }
    return resampledControlPoints;
}


const prevV = Vec3();
const tmpV1 = Vec3();
const tmpV2 = Vec3();
const tmpV3 = Vec3();

// easier to align to theses normals
function GetSmoothNormals(points: Vec3[]) {
    const nP: number = points.length;
    const smoothNormals: Vec3[] = [];
    if (points.length < 3) {
        for (let i = 0; i < points.length; ++i)
            smoothNormals.push(Vec3.normalize(Vec3(), points[i]));
        return smoothNormals;
    }
    let p0 = points[0];
    let p1 = points[1];
    let p2 = points[2];
    const p21 = Vec3.sub(tmpV1, p2, p1);
    const p01 =  Vec3.sub(tmpV2, p0, p1);
    const p0121 = Vec3.cross(tmpV3, p01, p21);
    Vec3.normalize(prevV, p0121);
    smoothNormals.push(Vec3.clone(prevV));
    for (let i = 1; i < points.length - 1; ++i) {
        p0 = points[i - 1];
        p1 = points[i];
        p2 = points[i + 1];
        const t = Vec3.normalize(tmpV1, Vec3.sub(tmpV1, p2, p0));
        const b = Vec3.normalize(tmpV2, Vec3.cross(tmpV2, t, prevV));
        const n = Vec3.normalize(Vec3(), Vec3.cross(tmpV3, t, b));
        Vec3.negate(n, n);
        Vec3.copy(prevV, n);
        smoothNormals.push(n);
    }
    const last = Vec3();
    Vec3.normalize(last, Vec3.cross(last,
        Vec3.sub(tmpV1, points[nP - 3], points[nP - 2]),
        Vec3.sub(tmpV2, points[nP - 2], points[nP - 1]))
    );
    smoothNormals.push(last);
    return smoothNormals;
}

const frameTmpV1 = Vec3();
const frameTmpV2 = Vec3();
const frameTmpV3 = Vec3();

function getFrame(reference: Vec3, tangent: Vec3) {
    const t = Vec3.normalize(Vec3(), tangent);
    // make reference vector orthogonal to tangent
    const proj_r_to_t = Vec3.scale(
        frameTmpV1, tangent, Vec3.dot(reference, tangent) / Vec3.dot(tangent, tangent)
    );
    const r = Vec3.normalize(Vec3(), Vec3.sub(frameTmpV2, reference, proj_r_to_t));
    // make bitangent vector orthogonal to the others
    const s = Vec3.normalize(Vec3(), Vec3.cross(frameTmpV3, t, r));
    return { t, r, s };
}

const mfTmpV1 = Vec3();
const mfTmpV2 = Vec3();
const mfTmpV3 = Vec3();
const mfTmpV4 = Vec3();
const mfTmpV5 = Vec3();
const mfTmpV6 = Vec3();
const mfTmpV7 = Vec3();
const mfTmpV8 = Vec3();
const mfTmpV9 = Vec3();

// easier to align to theses normals
// https://github.com/bzamecnik/gpg/blob/master/rotation-minimizing-frame/rmf.py
function GetMiniFrame(points: Vec3[], normals: Vec3[]) {
    const frames: Frame[] = [];
    const t0 = Vec3.normalize(mfTmpV1, Vec3.sub(mfTmpV1, points[1], points[0]));
    frames.push(getFrame(normals[0], t0));

    for (let i = 0; i < points.length - 2; ++i) {
        const t2 = Vec3.normalize(mfTmpV1, Vec3.sub(mfTmpV1, points[i + 2], points[i + 1]));
        const v1 = Vec3.sub(mfTmpV2, points[i + 1], points[i]); // this is tangeant
        const c1 = Vec3.dot(v1, v1);
        // compute r_i^L = R_1 * r_i
        const v1r = Vec3.scale(mfTmpV3, v1, (2.0 / c1) * Vec3.dot(v1, frames[i].r));
        const ref_L_i = Vec3.sub(mfTmpV4, frames[i].r, v1r);
        // compute t_i^L = R_1 * t_i
        const v1t = Vec3.scale(mfTmpV5, v1, (2.0 / c1) * Vec3.dot(v1, frames[i].t));
        const tan_L_i = Vec3.sub(mfTmpV6, frames[i].t, v1t);
        // # compute reflection vector of R_2
        const v2 =  Vec3.sub(mfTmpV7, t2, tan_L_i);
        const c2 = Vec3.dot(v2, v2);
        // compute r_(i+1) = R_2 * r_i^L
        const v2l = Vec3.scale(mfTmpV8, v1, (2.0 / c2) * Vec3.dot(v2, ref_L_i));
        const ref_next = Vec3.sub(mfTmpV9, ref_L_i, v2l); // ref_L_i - (2 / c2) * v2.dot(ref_L_i) * v2
        frames.push(getFrame(ref_next, t2)); // frames.append(Frame(ref_next, tangents[i+1]))
    }
    return frames;
}

const rpTmpVec1 = Vec3();

export function getMatFromResamplePoints(points: NumberArray, segmentLength: number, resample: boolean) {
    let new_points: Vec3[] = [];
    if (resample) new_points = ResampleControlPoints(points, segmentLength);
    else {
        for (let idx = 0; idx < points.length / 3; ++idx){
            new_points.push(Vec3.fromArray(Vec3.zero(), points, idx * 3));
        }
    }
    const npoints = new_points.length;
    const new_normal = GetSmoothNormals(new_points);
    const frames = GetMiniFrame(new_points, new_normal);
    const limit = npoints;
    const transforms: Mat4[] = [];
    const pti = Vec3.copy(rpTmpVec1, new_points[0]);
    for (let i = 0; i < npoints - 2; ++i) {
        const pti1: Vec3 = new_points[i + 1]; // Vec3.create(points[(i+1)*3],points[(i+1)*3+1],points[(i+1)*3+2]);
        const d = Vec3.distance(pti, pti1);
        if (d >= segmentLength) {
            // use twist or random?
            const quat = Quat.rotationTo(Quat.zero(), Vec3.create(0, 0, 1), frames[i].t); // Quat.rotationTo(Quat.zero(), Vec3.create(0,0,1),new_normal[i]);//Quat.rotationTo(Quat.zero(), Vec3.create(0,0,1),direction);new_normal
            const rq = Quat.setAxisAngle(Quat.zero(), frames[i].t, Math.random() * 3.60 ); // Quat.setAxisAngle(Quat.zero(),direction, Math.random()*3.60 );//Quat.identity();//
            const m = Mat4.fromQuat(Mat4.zero(), Quat.multiply(Quat.zero(), rq, quat)); // Mat4.fromQuat(Mat4.zero(),Quat.multiply(Quat.zero(),quat1,quat2));//Mat4.fromQuat(Mat4.zero(),quat);//Mat4.identity();//Mat4.fromQuat(Mat4.zero(),Quat.multiply(Quat.zero(),rq,quat));
            // let pos:Vec3 = Vec3.add(Vec3.zero(),pti1,pti)
            // pos = Vec3.scale(pos,pos,1.0/2.0);
            // Vec3.makeRotation(Mat4.zero(),Vec3.create(0,0,1),frames[i].t);//
            Mat4.setTranslation(m, pti1);
            // let m2:Mat4 = GetTubePropertiesMatrix(pti,pti1);
            // let q:Quat = Quat.rotationTo(Quat.zero(), Vec3.create(0,1,0),Vec3.create(0,0,1))
            // m2=Mat4.mul(Mat4.identity(),Mat4.fromQuat(Mat4.zero(),q),m2);
            transforms.push(m);
            Vec3.copy(pti, pti1);
        }
        if (transforms.length >= limit) break;
    }
    return transforms;
}