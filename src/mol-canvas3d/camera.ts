/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3, Vec4, EPSILON } from 'mol-math/linear-algebra'
import { Viewport, cameraProject, cameraUnproject } from './camera/util';
import { Object3D } from 'mol-gl/object3d';
import { BehaviorSubject } from 'rxjs';
import { CameraTransitionManager } from './camera/transition';

export { Camera }

// TODO: slab controls that modify near/far planes?

class Camera implements Object3D {
    readonly updatedViewProjection = new BehaviorSubject<Camera>(this);

    readonly view: Mat4 = Mat4.identity();
    readonly projection: Mat4 = Mat4.identity();
    readonly projectionView: Mat4 = Mat4.identity();
    readonly inverseProjectionView: Mat4 = Mat4.identity();

    readonly viewport: Viewport;
    readonly state: Readonly<Camera.Snapshot> = Camera.createDefaultSnapshot();
    readonly viewOffset: Camera.ViewOffset = {
        enabled: false,
        fullWidth: 1, fullHeight: 1,
        offsetX: 0, offsetY: 0,
        width: 1, height: 1
    }

    readonly transition: CameraTransitionManager = new CameraTransitionManager(this);

    get position() { return this.state.position; }
    set position(v: Vec3) { Vec3.copy(this.state.position, v); }

    get direction() { return this.state.direction; }
    set direction(v: Vec3) { Vec3.copy(this.state.direction, v); }

    get up() { return this.state.up; }
    set up(v: Vec3) { Vec3.copy(this.state.up, v); }

    get target() { return this.state.target; }
    set target(v: Vec3) { Vec3.copy(this.state.target, v); }

    private prevProjection = Mat4.identity();
    private prevView = Mat4.identity();
    private deltaDirection = Vec3.zero();
    private newPosition = Vec3.zero();

    updateMatrices() {
        const snapshot = this.state as Camera.Snapshot;
        const height = 2 * Math.tan(snapshot.fov / 2) * Vec3.distance(snapshot.position, snapshot.target);
        snapshot.zoom = this.viewport.height / height;

        switch (this.state.mode) {
            case 'orthographic': updateOrtho(this); break;
            case 'perspective': updatePers(this); break;
            default: throw new Error('unknown camera mode');
        }

        const changed = !Mat4.areEqual(this.projection, this.prevProjection, EPSILON.Value) || !Mat4.areEqual(this.view, this.prevView, EPSILON.Value);

        Mat4.mul(this.projectionView, this.projection, this.view)
        Mat4.invert(this.inverseProjectionView, this.projectionView)


        if (changed) {
            Mat4.mul(this.projectionView, this.projection, this.view)
            Mat4.invert(this.inverseProjectionView, this.projectionView)

            Mat4.copy(this.prevView, this.view);
            Mat4.copy(this.prevProjection, this.projection);
            this.updatedViewProjection.next(this);
        }

        return changed;
    }

    setState(snapshot: Partial<Camera.Snapshot>) {
        this.transition.apply(snapshot);
    }

    getSnapshot() {
        const ret = Camera.createDefaultSnapshot();
        Camera.copySnapshot(ret, this.state);
        return ret;
    }

    getFocus(target: Vec3, radius: number): Partial<Camera.Snapshot> {
        const fov = this.state.fov
        const { width, height } = this.viewport
        const aspect = width / height
        const aspectFactor = (height < width ? 1 : aspect)
        const currentDistance = Vec3.distance(this.state.position, target)
        const targetDistance = Math.abs((radius / aspectFactor) / Math.sin(fov / 2))
        const deltaDistance = Math.abs(currentDistance - targetDistance)

        Vec3.sub(this.deltaDirection, this.state.position, target)
        Vec3.setMagnitude(this.deltaDirection, this.state.direction, deltaDistance)
        if (currentDistance < targetDistance) Vec3.negate(this.deltaDirection, this.deltaDirection)
        Vec3.add(this.newPosition, this.state.position, this.deltaDirection)

        return { target, position: Vec3.clone(this.newPosition) };
    }

    focus(target: Vec3, radius: number) {
        this.setState(this.getFocus(target, radius));
    }

    // lookAt(target: Vec3) {
    //     cameraLookAt(this.position, this.up, this.direction, target);
    // }

    // translate(v: Vec3) {
    //     Vec3.add(this.position, this.position, v);
    //     cameraLookAt(this.position, this.up, this.direction, this.target);
    // }

    project(out: Vec4, point: Vec3) {
        return cameraProject(out, point, this.viewport, this.projectionView)
    }

    unproject(out: Vec3, point: Vec3) {
        return cameraUnproject(out, point, this.viewport, this.inverseProjectionView)
    }

    dispose() {
        this.updatedViewProjection.complete();
    }

    constructor(state?: Partial<Camera.Snapshot>, viewport = Viewport.create(-1, -1, 1, 1)) {
        this.viewport = viewport;
        Camera.copySnapshot(this.state, state);
    }

}

namespace Camera {
    export type Mode = 'perspective' | 'orthographic'

    export interface ClippingInfo {
        near: number,
        far: number,
        fogNear: number,
        fogFar: number
    }

    /**
     * Sets an offseted view in a larger frustum. This is useful for
     * - multi-window or multi-monitor/multi-machine setups
     * - jittering the camera position for
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

    export function setViewOffset(out: ViewOffset, fullWidth: number, fullHeight: number, offsetX: number, offsetY: number, width: number, height: number) {
        out.fullWidth = fullWidth
        out.fullHeight = fullHeight
        out.offsetX = offsetX
        out.offsetY = offsetY
        out.width = width
        out.height = height
    }

    export function createDefaultSnapshot(): Snapshot {
        return {
            mode: 'perspective',

            position: Vec3.zero(),
            direction: Vec3.create(0, 0, -1),
            up: Vec3.create(0, 1, 0),

            target: Vec3.create(0, 0, 0),

            near: 1,
            far: 10000,
            fogNear: 1,
            fogFar: 10000,

            fov: Math.PI / 4,
            zoom: 1,
        };
    }

    export interface Snapshot {
        mode: Mode,

        position: Vec3,
        // Normalized camera direction
        direction: Vec3,
        up: Vec3,
        target: Vec3,

        near: number,
        far: number,
        fogNear: number,
        fogFar: number,

        fov: number,
        zoom: number,
    }

    export function copySnapshot(out: Snapshot, source?: Partial<Snapshot>) {
        if (!source) return;

        if (typeof source.mode !== 'undefined') out.mode = source.mode;

        if (typeof source.position !== 'undefined') Vec3.copy(out.position, source.position);
        if (typeof source.direction !== 'undefined') Vec3.copy(out.direction, source.direction);
        if (typeof source.up !== 'undefined') Vec3.copy(out.up, source.up);
        if (typeof source.target !== 'undefined') Vec3.copy(out.target, source.target);

        if (typeof source.near !== 'undefined') out.near = source.near;
        if (typeof source.far !== 'undefined') out.far = source.far;
        if (typeof source.fogNear !== 'undefined') out.fogNear = source.fogNear;
        if (typeof source.fogFar !== 'undefined') out.fogFar = source.fogFar;

        if (typeof source.fov !== 'undefined') out.fov = source.fov;
        if (typeof source.zoom !== 'undefined') out.zoom = source.zoom;
    }
}

const _center = Vec3.zero();
function updateOrtho(camera: Camera) {
    const { viewport, state: { zoom, near, far }, viewOffset } = camera

    const fullLeft = -(viewport.width - viewport.x) / 2
    const fullRight = (viewport.width - viewport.x) / 2
    const fullTop = (viewport.height - viewport.y) / 2
    const fullBottom = -(viewport.height - viewport.y) / 2

    const dx = (fullRight - fullLeft) / (2 * zoom)
    const dy = (fullTop - fullBottom) / (2 * zoom)
    const cx = (fullRight + fullLeft) / 2
    const cy = (fullTop + fullBottom) / 2

    let left = cx - dx
    let right = cx + dx
    let top = cy + dy
    let bottom = cy - dy

    if (viewOffset && viewOffset.enabled) {
        const zoomW = zoom / (viewOffset.width / viewOffset.fullWidth)
        const zoomH = zoom / (viewOffset.height / viewOffset.fullHeight)
        const scaleW = (fullRight - fullLeft) / viewOffset.width
        const scaleH = (fullTop - fullBottom) / viewOffset.height
        left += scaleW * (viewOffset.offsetX / zoomW)
        right = left + scaleW * (viewOffset.width / zoomW)
        top -= scaleH * (viewOffset.offsetY / zoomH)
        bottom = top - scaleH * (viewOffset.height / zoomH)
    }

    // build projection matrix
    Mat4.ortho(camera.projection, left, right, top, bottom, near, far)

    // build view matrix
    Vec3.add(_center, camera.position, camera.direction)
    Mat4.lookAt(camera.view, camera.position, _center, camera.up)
}

function updatePers(camera: Camera) {
    const aspect = camera.viewport.width / camera.viewport.height

    const { state: { fov, near, far }, viewOffset } = camera

    let top = near * Math.tan(0.5 * fov)
    let height = 2 * top
    let width = aspect * height
    let left = -0.5 * width

    if (viewOffset && viewOffset.enabled) {
        left += viewOffset.offsetX * width / viewOffset.fullWidth
        top -= viewOffset.offsetY * height / viewOffset.fullHeight
        width *= viewOffset.width / viewOffset.fullWidth
        height *= viewOffset.height / viewOffset.fullHeight
    }

    // build projection matrix
    Mat4.perspective(camera.projection, left, left + width, top, top - height, near, far)

    // build view matrix
    Vec3.add(_center, camera.position, camera.direction)
    Mat4.lookAt(camera.view, camera.position, _center, camera.up)
}