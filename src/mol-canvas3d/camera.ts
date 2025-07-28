/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3, Vec4, EPSILON } from '../mol-math/linear-algebra';
import { Viewport, cameraProject, cameraUnproject } from './camera/util';
import { CameraTransitionManager } from './camera/transition';
import { BehaviorSubject } from 'rxjs';
import { Scene } from '../mol-gl/scene';
import { assertUnreachable } from '../mol-util/type-helpers';
import { Ray3D } from '../mol-math/geometry/primitives/ray3d';

export type { ICamera };

interface ICamera {
    readonly viewport: Viewport,
    readonly view: Mat4,
    readonly projection: Mat4,
    readonly projectionView: Mat4,
    readonly inverseProjectionView: Mat4,
    readonly state: Readonly<Camera.Snapshot>,
    readonly viewOffset: Camera.ViewOffset,
    readonly far: number,
    readonly near: number,
    readonly fogFar: number,
    readonly fogNear: number,
    readonly headRotation: Mat4,
}

const tmpClip = Vec4();

export class Camera implements ICamera {
    readonly view: Mat4 = Mat4.identity();
    readonly projection: Mat4 = Mat4.identity();
    readonly projectionView: Mat4 = Mat4.identity();
    readonly inverseProjectionView: Mat4 = Mat4.identity();
    readonly headRotation: Mat4 = Mat4.zero();

    readonly viewport: Viewport;
    readonly state: Readonly<Camera.Snapshot> = Camera.createDefaultSnapshot();
    readonly viewOffset = Camera.ViewOffset();

    near = 1;
    far = 10000;
    fogNear = 5000;
    fogFar = 10000;
    zoom = 1;

    readonly transition: CameraTransitionManager = new CameraTransitionManager(this);
    readonly stateChanged = new BehaviorSubject<Partial<Camera.Snapshot>>(this.state);

    get position() { return this.state.position; }
    set position(v: Vec3) { Vec3.copy(this.state.position, v); }

    get up() { return this.state.up; }
    set up(v: Vec3) { Vec3.copy(this.state.up, v); }

    get target() { return this.state.target; }
    set target(v: Vec3) { Vec3.copy(this.state.target, v); }

    private prevProjection = Mat4.identity();
    private prevView = Mat4.identity();
    private deltaDirection = Vec3();
    private newPosition = Vec3();

    update() {
        const snapshot = this.state as Camera.Snapshot;
        if (snapshot.radiusMax === 0) {
            return false;
        }

        const height = 2 * Math.tan(snapshot.fov / 2) * Vec3.distance(snapshot.position, snapshot.target) * this.state.scale;
        this.zoom = this.viewport.height / height;

        updateClip(this);

        switch (this.state.mode) {
            case 'orthographic': updateOrtho(this); break;
            case 'perspective': updatePers(this); break;
            default: assertUnreachable(this.state.mode);
        }

        const changed = !Mat4.areEqual(this.projection, this.prevProjection, EPSILON) || !Mat4.areEqual(this.view, this.prevView, EPSILON);

        if (changed) {
            Mat4.mul(this.projectionView, this.projection, this.view);
            if (!Mat4.tryInvert(this.inverseProjectionView, this.projectionView)) {
                Mat4.copy(this.view, this.prevView);
                Mat4.copy(this.projection, this.prevProjection);
                Mat4.mul(this.projectionView, this.projection, this.view);
                return false;
            }

            Mat4.copy(this.prevView, this.view);
            Mat4.copy(this.prevProjection, this.projection);
        }

        return changed;
    }

    setState(snapshot: Partial<Camera.Snapshot>, durationMs?: number) {
        this.transition.apply(snapshot, durationMs);
        this.stateChanged.next(snapshot);
    }

    getSnapshot() {
        return Camera.copySnapshot(Camera.createDefaultSnapshot(), this.state);
    }

    getTargetDistance(radius: number) {
        return Camera.targetDistance(radius, this.state.mode, this.state.fov, this.viewport.width, this.viewport.height);
    }

    getFocus(target: Vec3, radius: number, up?: Vec3, dir?: Vec3, snapshot?: Partial<Camera.Snapshot>): Partial<Camera.Snapshot> {
        const r = Math.max(radius, 0.01);
        const targetDistance = this.getTargetDistance(r);

        Vec3.sub(this.deltaDirection, snapshot?.target ?? this.target, snapshot?.position ?? this.position);
        if (dir) Vec3.matchDirection(this.deltaDirection, dir, this.deltaDirection);
        Vec3.setMagnitude(this.deltaDirection, this.deltaDirection, targetDistance);
        Vec3.sub(this.newPosition, target, this.deltaDirection);

        const state = Camera.copySnapshot(Camera.createDefaultSnapshot(), this.state);
        state.target = Vec3.clone(target);
        state.radius = r;
        state.position = Vec3.clone(this.newPosition);
        if (up) Vec3.matchDirection(state.up, up, state.up);

        return state;
    }

    getCenter(target: Vec3, radius?: number): Partial<Camera.Snapshot> {
        Vec3.sub(this.deltaDirection, this.target, this.position);
        Vec3.sub(this.newPosition, target, this.deltaDirection);

        const state = Camera.copySnapshot(Camera.createDefaultSnapshot(), this.state);
        state.target = Vec3.clone(target);
        state.position = Vec3.clone(this.newPosition);
        if (radius) state.radius = Math.max(radius, 0.01);

        return state;
    }

    getInvariantFocus(target: Vec3, radius: number, up: Vec3, dir: Vec3): Partial<Camera.Snapshot> {
        const r = Math.max(radius, 0.01);
        const targetDistance = this.getTargetDistance(r);

        Vec3.copy(this.deltaDirection, dir);
        Vec3.setMagnitude(this.deltaDirection, this.deltaDirection, targetDistance);
        Vec3.sub(this.newPosition, target, this.deltaDirection);

        const state = Camera.copySnapshot(Camera.createDefaultSnapshot(), this.state);
        state.target = Vec3.clone(target);
        state.radius = r;
        state.position = Vec3.clone(this.newPosition);
        Vec3.copy(state.up, up);

        return state;
    }

    focus(target: Vec3, radius: number, durationMs?: number, up?: Vec3, dir?: Vec3) {
        if (radius > 0) {
            this.setState(this.getFocus(target, radius, up, dir), durationMs);
        }
    }

    center(target: Vec3, durationMs?: number) {
        this.setState(this.getCenter(target), durationMs);
    }

    /** Transform point into 2D window coordinates. */
    project(out: Vec4, point: Vec3) {
        return cameraProject(out, point, this.viewport, this.projectionView);
    }

    /**
     * Transform point from screen space to 3D coordinates.
     * The point must have `x` and `y` set to 2D window coordinates
     * and `z` between 0 (near) and 1 (far); the optional `w` is not used.
     */
    unproject(out: Vec3, point: Vec3 | Vec4) {
        return cameraUnproject(out, point, this.viewport, this.inverseProjectionView);
    }

    /** World space pixel size at given `point` */
    getPixelSize(point: Vec3) {
        this.project(tmpClip, point);
        const w = tmpClip[3];
        const rx = this.viewport.width;
        const P00 = this.projection[0];
        return (2 / w) / (rx * Math.abs(P00));
    }

    getRay(out: Ray3D, x: number, y: number) {
        if (this.state.mode === 'orthographic') {
            Vec3.set(out.origin, x, y, 0);
            this.unproject(out.origin, out.origin);
            Vec3.normalize(out.direction, Vec3.sub(out.direction, this.target, this.position));
            Vec3.scaleAndAdd(out.origin, out.origin, out.direction, -this.near);
        } else {
            Vec3.copy(out.origin, this.state.position);
            Vec3.scale(out.origin, out.origin, this.state.scale);
            Vec3.set(out.direction, x, y, 0.5);
            this.unproject(out.direction, out.direction);
            Vec3.normalize(out.direction, Vec3.sub(out.direction, out.direction, out.origin));
        }
        return out;
    }

    constructor(state?: Partial<Camera.Snapshot>, viewport = Viewport.create(0, 0, 128, 128)) {
        this.viewport = viewport;
        Camera.copySnapshot(this.state, state);
    }
}

export namespace Camera {
    export type Mode = 'perspective' | 'orthographic'

    export type SnapshotProvider = Partial<Snapshot> | ((scene: Scene, camera: Camera) => Partial<Snapshot>)

    /**
     * Sets an offseted view in a larger frustum. This is useful for
     * - multi-window or multi-monitor/multi-machine setups
     * - jittering the camera position for sampling
     */
    export interface ViewOffset {
        enabled: boolean,
        fullWidth: number,
        fullHeight: number,
        offsetX: number,
        offsetY: number,
        width: number,
        height: number
    }

    export function ViewOffset(): ViewOffset {
        return {
            enabled: false,
            fullWidth: 1, fullHeight: 1,
            offsetX: 0, offsetY: 0,
            width: 1, height: 1
        };
    }

    export function setViewOffset(out: ViewOffset, fullWidth: number, fullHeight: number, offsetX: number, offsetY: number, width: number, height: number) {
        out.fullWidth = fullWidth;
        out.fullHeight = fullHeight;
        out.offsetX = offsetX;
        out.offsetY = offsetY;
        out.width = width;
        out.height = height;
    }

    export function copyViewOffset(out: ViewOffset, view: ViewOffset) {
        out.enabled = view.enabled;
        out.fullWidth = view.fullWidth;
        out.fullHeight = view.fullHeight;
        out.offsetX = view.offsetX;
        out.offsetY = view.offsetY;
        out.width = view.width;
        out.height = view.height;
    }

    export function targetDistance(radius: number, mode: Mode, fov: number, width: number, height: number) {
        const r = Math.max(radius, 0.01);
        const aspect = width / height;
        const aspectFactor = (height < width ? 1 : aspect);
        if (mode === 'orthographic')
            return Math.abs((r / aspectFactor) / Math.tan(fov / 2));
        else
            return Math.abs((r / aspectFactor) / Math.sin(fov / 2));
    }

    export function createDefaultSnapshot(): Snapshot {
        return {
            mode: 'perspective',
            fov: Math.PI / 4,

            position: Vec3.create(0, 0, 100),
            up: Vec3.create(0, 1, 0),
            target: Vec3.create(0, 0, 0),

            radius: 0,
            radiusMax: 10,
            fog: 50,
            clipFar: true,
            minNear: 5,
            minFar: 0,

            scale: 1,
        };
    }

    export interface Snapshot {
        mode: Mode
        fov: number

        position: Vec3
        up: Vec3
        target: Vec3

        radius: number
        radiusMax: number
        fog: number
        clipFar: boolean
        minNear: number
        minFar: number

        scale: number
    }

    export function copySnapshot(out: Snapshot, source?: Partial<Snapshot>) {
        if (!source) return out;

        if (typeof source.mode !== 'undefined') out.mode = source.mode;
        if (typeof source.fov !== 'undefined') out.fov = source.fov;

        if (typeof source.position !== 'undefined') Vec3.copy(out.position, source.position);
        if (typeof source.up !== 'undefined') Vec3.copy(out.up, source.up);
        if (typeof source.target !== 'undefined') Vec3.copy(out.target, source.target);

        if (typeof source.radius !== 'undefined') out.radius = source.radius;
        if (typeof source.radiusMax !== 'undefined') out.radiusMax = source.radiusMax;
        if (typeof source.fog !== 'undefined') out.fog = source.fog;
        if (typeof source.clipFar !== 'undefined') out.clipFar = source.clipFar;
        if (typeof source.minNear !== 'undefined') out.minNear = source.minNear;
        if (typeof source.minFar !== 'undefined') out.minFar = source.minFar;

        if (typeof source.scale !== 'undefined') out.scale = source.scale;

        return out;
    }

    export function areSnapshotsEqual(a: Snapshot, b: Snapshot) {
        return a.mode === b.mode
            && a.fov === b.fov
            && a.radius === b.radius
            && a.radiusMax === b.radiusMax
            && a.fog === b.fog
            && a.clipFar === b.clipFar
            && a.minNear === b.minNear
            && a.minFar === b.minFar
            && a.scale === b.scale
            && Vec3.exactEquals(a.position, b.position)
            && Vec3.exactEquals(a.up, b.up)
            && Vec3.exactEquals(a.target, b.target);
    }
}

const tmpPosition = Vec3();
const tmpTarget = Vec3();

function updateView(camera: Camera) {
    if (camera.state.scale === 1) {
        Mat4.lookAt(camera.view, camera.state.position, camera.state.target, camera.state.up);
    } else {
        Vec3.scale(tmpPosition, camera.state.position, camera.state.scale);
        Vec3.scale(tmpTarget, camera.state.target, camera.state.scale);
        Mat4.lookAt(camera.view, tmpPosition, tmpTarget, camera.state.up);
    }
}

function updateOrtho(camera: Camera) {
    const { viewport, zoom, near, far, viewOffset } = camera;

    const fullLeft = -viewport.width / 2;
    const fullRight = viewport.width / 2;
    const fullTop = viewport.height / 2;
    const fullBottom = -viewport.height / 2;

    const dx = (fullRight - fullLeft) / (2 * zoom);
    const dy = (fullTop - fullBottom) / (2 * zoom);
    const cx = (fullRight + fullLeft) / 2;
    const cy = (fullTop + fullBottom) / 2;

    let left = cx - dx;
    let right = cx + dx;
    let top = cy + dy;
    let bottom = cy - dy;

    if (viewOffset.enabled) {
        const zoomW = zoom / (viewOffset.width / viewOffset.fullWidth);
        const zoomH = zoom / (viewOffset.height / viewOffset.fullHeight);
        const scaleW = (fullRight - fullLeft) / viewOffset.width;
        const scaleH = (fullTop - fullBottom) / viewOffset.height;
        left += scaleW * (viewOffset.offsetX / zoomW);
        right = left + scaleW * (viewOffset.width / zoomW);
        top -= scaleH * (viewOffset.offsetY / zoomH);
        bottom = top - scaleH * (viewOffset.height / zoomH);
    }

    // build projection matrix
    Mat4.ortho(camera.projection, left, right, top, bottom, near, far);

    // build view matrix
    updateView(camera);
}

function updatePers(camera: Camera) {
    const aspect = camera.viewport.width / camera.viewport.height;

    const { near, far, viewOffset } = camera;

    let top = near * Math.tan(0.5 * camera.state.fov);
    let height = 2 * top;
    let width = aspect * height;
    let left = -0.5 * width;

    if (viewOffset.enabled) {
        left += viewOffset.offsetX * width / viewOffset.fullWidth;
        top -= viewOffset.offsetY * height / viewOffset.fullHeight;
        width *= viewOffset.width / viewOffset.fullWidth;
        height *= viewOffset.height / viewOffset.fullHeight;
    }

    // build projection matrix
    Mat4.perspective(camera.projection, left, left + width, top, top - height, near, far);

    // build view matrix
    updateView(camera);
}

function updateClip(camera: Camera) {
    let { radius, radiusMax, mode, fog, clipFar, minNear, minFar, scale } = camera.state;
    radiusMax *= scale;
    minFar *= scale;
    minNear *= scale;
    radius *= scale;

    const minRadius = 0.01 * scale;
    if (radius < minRadius) radius = minRadius;

    const normalizedFar = Math.max(clipFar ? radius : radiusMax, minFar);
    Vec3.scale(tmpTarget, camera.state.target, scale);
    Vec3.scale(tmpPosition, camera.state.position, scale);
    const cameraDistance = Vec3.distance(tmpPosition, tmpTarget);
    let near = cameraDistance - radius;
    let far = cameraDistance + normalizedFar;

    if (mode === 'perspective') {
        // set at least to 5 to avoid slow sphere impostor rendering
        near = Math.max(Math.min(radiusMax, minNear), near);
        far = Math.max(minNear, far);
    } else {
        // not too close to 0 as it causes issues with outline rendering
        near = Math.max(Math.min(radiusMax, minNear), near);
        far = Math.max(minNear, far);
    }

    if (near === far) {
        // make sure near and far are not identical to avoid Infinity in the projection matrix
        far = near + 0.01 * scale;
    }

    const fogNearFactor = -(50 - fog) / 50;
    const fogNear = cameraDistance - (normalizedFar * fogNearFactor);
    const fogFar = far;

    camera.near = near;
    camera.far = far;
    camera.fogNear = fogNear;
    camera.fogFar = fogFar;
}