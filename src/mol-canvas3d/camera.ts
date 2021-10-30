/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3, Vec4, EPSILON } from '../mol-math/linear-algebra';
import { Viewport, cameraProject, cameraUnproject } from './camera/util';
import { CameraTransitionManager } from './camera/transition';
import { BehaviorSubject } from 'rxjs';
import { Scene } from '../mol-gl/scene';

export { ICamera, Camera };

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
}

const tmpPos1 = Vec3();
const tmpPos2 = Vec3();
const tmpClip = Vec4();

class Camera implements ICamera {
    readonly view: Mat4 = Mat4.identity();
    readonly projection: Mat4 = Mat4.identity();
    readonly projectionView: Mat4 = Mat4.identity();
    readonly inverseProjectionView: Mat4 = Mat4.identity();

    private pixelScale: number
    get pixelRatio() {
        const dpr = (typeof window !== 'undefined') ? window.devicePixelRatio : 1;
        return dpr * this.pixelScale;
    }

    readonly viewport: Viewport;
    readonly state: Readonly<Camera.Snapshot> = Camera.createDefaultSnapshot();
    readonly viewOffset = Camera.ViewOffset();

    near = 1
    far = 10000
    fogNear = 5000
    fogFar = 10000
    zoom = 1

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

        const height = 2 * Math.tan(snapshot.fov / 2) * Vec3.distance(snapshot.position, snapshot.target);
        this.zoom = this.viewport.height / height;

        updateClip(this);

        switch (this.state.mode) {
            case 'orthographic': updateOrtho(this); break;
            case 'perspective': updatePers(this); break;
            default: throw new Error('unknown camera mode');
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
        return Camera.targetDistance(radius, this.state.fov, this.viewport.width, this.viewport.height);
    }

    getFocus(target: Vec3, radius: number, up?: Vec3, dir?: Vec3): Partial<Camera.Snapshot> {
        const r = Math.max(radius, 0.01);
        const targetDistance = this.getTargetDistance(r);

        Vec3.sub(this.deltaDirection, this.target, this.position);
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
        // project -> unproject of `point` does not exactly return the same
        // to get a sufficiently accurate measure we unproject the original
        // clip position in addition to the one shifted bey one pixel
        this.project(tmpClip, point);
        this.unproject(tmpPos1, tmpClip);
        tmpClip[0] += 1;
        this.unproject(tmpPos2, tmpClip);
        return Vec3.distance(tmpPos1, tmpPos2);
    }

    constructor(state?: Partial<Camera.Snapshot>, viewport = Viewport.create(0, 0, 128, 128), props: Partial<{ pixelScale: number }> = {}) {
        this.viewport = viewport;
        this.pixelScale = props.pixelScale || 1;
        Camera.copySnapshot(this.state, state);
    }
}

namespace Camera {
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

    export function targetDistance(radius: number, fov: number, width: number, height: number) {
        const r = Math.max(radius, 0.01);
        const aspect = width / height;
        const aspectFactor = (height < width ? 1 : aspect);
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
            clipFar: true
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

        return out;
    }

    export function areSnapshotsEqual(a: Snapshot, b: Snapshot) {
        return a.mode === b.mode
            && a.fov === b.fov
            && a.radius === b.radius
            && a.radiusMax === b.radiusMax
            && a.fog === b.fog
            && a.clipFar === b.clipFar
            && Vec3.exactEquals(a.position, b.position)
            && Vec3.exactEquals(a.up, b.up)
            && Vec3.exactEquals(a.target, b.target);
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
    Mat4.lookAt(camera.view, camera.position, camera.target, camera.up);
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
    Mat4.lookAt(camera.view, camera.position, camera.target, camera.up);
}

function updateClip(camera: Camera) {
    let { radius, radiusMax, mode, fog, clipFar } = camera.state;
    if (radius < 0.01) radius = 0.01;

    const normalizedFar = clipFar ? radius : radiusMax;
    const cameraDistance = Vec3.distance(camera.position, camera.target);
    let near = cameraDistance - radius;
    let far = cameraDistance + normalizedFar;

    const fogNearFactor = -(50 - fog) / 50;
    const fogNear = cameraDistance - (normalizedFar * fogNearFactor);
    const fogFar = far;

    if (mode === 'perspective') {
        // set at least to 5 to avoid slow sphere impostor rendering
        near = Math.max(Math.min(radiusMax, 5), near);
        far = Math.max(5, far);
    } else {
        // not too close to 0 as it causes issues with outline rendering
        near = Math.max(Math.min(radiusMax, 5), near);
        far = Math.max(5, far);
    }

    if (near === far) {
        // make sure near and far are not identical to avoid Infinity in the projection matrix
        far = near + 0.01;
    }

    camera.near = near;
    camera.far = 2 * far; // avoid precision issues distingushing far objects from background
    camera.fogNear = fogNear;
    camera.fogFar = fogFar;
}