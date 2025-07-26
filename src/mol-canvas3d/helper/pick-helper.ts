/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderer } from '../../mol-gl/renderer';
import { Scene } from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { spiral2d } from '../../mol-math/misc';
import { isTimingMode } from '../../mol-util/debug';
import { Camera } from '../camera';
import { StereoCamera } from '../camera/stereo';
import { cameraUnproject, Viewport } from '../camera/util';
import { Helper } from '../helper/helper';
import { AsyncPickData, AsyncPickStatus, checkAsyncPickingSupport, PickBuffers, PickData, PickOptions, PickPass } from '../passes/pick';

export class PickHelper {
    dirty = true;

    private pickPadding: number;
    private buffers = new PickBuffers(this.webgl, this.pickPass);
    private viewport = Viewport();

    private pickRatio: number;
    private pickX: number;
    private pickY: number;
    private pickWidth: number;
    private pickHeight: number;
    private halfPickWidth: number;

    private spiral: [number, number][];

    setViewport(x: number, y: number, width: number, height: number) {
        Viewport.set(this.viewport, x, y, width, height);
        this.update();
    }

    setPickPadding(pickPadding: number) {
        if (this.pickPadding !== pickPadding) {
            this.pickPadding = pickPadding;
            this.update();
        }
    }

    private update() {
        const { x, y, width, height } = this.viewport;

        this.pickRatio = this.pickPass.pickRatio;
        this.pickX = Math.ceil(x * this.pickRatio);
        this.pickY = Math.ceil(y * this.pickRatio);

        const pickWidth = Math.floor(width * this.pickRatio);
        const pickHeight = Math.floor(height * this.pickRatio);

        if (pickWidth !== this.pickWidth || pickHeight !== this.pickHeight) {
            this.pickWidth = pickWidth;
            this.pickHeight = pickHeight;
            this.halfPickWidth = Math.floor(this.pickWidth / 2);

            this.buffers.setViewport(this.pickX, this.pickY, this.pickWidth, this.pickHeight);
        }

        this.spiral = spiral2d(Math.ceil(this.pickRatio * this.pickPadding));
        this.dirty = true;
    }

    private render(camera: Camera | StereoCamera) {
        if (isTimingMode) this.webgl.timer.mark('PickHelper.render', { captureStats: true });
        const { pickX, pickY, pickWidth, pickHeight, halfPickWidth } = this;
        const { renderer, scene, helper } = this;

        renderer.setTransparentBackground(false);
        renderer.setDrawingBufferSize(pickWidth, pickHeight);
        renderer.setPixelRatio(this.pickRatio);

        if (StereoCamera.is(camera)) {
            renderer.setViewport(pickX, pickY, halfPickWidth, pickHeight);
            this.pickPass.render(renderer, camera.left, scene, helper);

            renderer.setViewport(pickX + halfPickWidth, pickY, pickWidth - halfPickWidth, pickHeight);
            this.pickPass.render(renderer, camera.right, scene, helper);
        } else {
            renderer.setViewport(pickX, pickY, pickWidth, pickHeight);
            this.pickPass.render(renderer, camera, scene, helper);
        }

        this.dirty = false;
        if (isTimingMode) this.webgl.timer.markEnd('PickHelper.render');
    }

    private identifyInternal(x: number, y: number, camera: Camera | StereoCamera): PickData | undefined {
        if (this.webgl.isContextLost) return;

        const { webgl, pickRatio } = this;
        if (webgl.isContextLost) return;

        x *= webgl.pixelRatio;
        y *= webgl.pixelRatio;
        y = this.pickPass.drawingBufferHeight - y; // flip y

        const { viewport } = this;

        // check if within viewport
        if (x < viewport.x ||
            y < viewport.y ||
            x > viewport.x + viewport.width ||
            y > viewport.y + viewport.height
        ) return;

        const xv = x - viewport.x;
        const yv = y - viewport.y;

        const xp = Math.floor(xv * pickRatio);
        const yp = Math.floor(yv * pickRatio);

        const pickingId = this.buffers.getPickingId(xp, yp);
        if (pickingId === undefined) return;

        const z = this.buffers.getDepth(xp, yp);
        const position = Vec3.create(x, y, z);
        if (StereoCamera.is(camera)) {
            const halfWidth = Math.floor(viewport.width / 2);
            if (x > viewport.x + halfWidth) {
                position[0] = viewport.x + (xv - halfWidth) * 2;
                cameraUnproject(position, position, viewport, camera.right.inverseProjectionView);
            } else {
                position[0] = viewport.x + xv * 2;
                cameraUnproject(position, position, viewport, camera.left.inverseProjectionView);
            }
        } else {
            cameraUnproject(position, position, viewport, camera.inverseProjectionView);
        }

        return { id: pickingId, position };
    }

    private prepare() {
        if (this.pickRatio !== this.pickPass.pickRatio) {
            this.update();
        }
    }

    private getPickData(x: number, y: number, camera: Camera | StereoCamera): PickData | undefined {
        for (const d of this.spiral) {
            const pickData = this.identifyInternal(x + d[0], y + d[1], camera);
            if (pickData) return pickData;
        }
    }

    identify(x: number, y: number, camera: Camera | StereoCamera): PickData | undefined {
        this.prepare();

        if (this.dirty) {
            if (isTimingMode) this.webgl.timer.mark('PickHelper.identify');
            this.render(camera);
            this.buffers.read();
            if (isTimingMode) this.webgl.timer.markEnd('PickHelper.identify');
        }

        return this.getPickData(x, y, camera);
    }

    asyncIdentify(x: number, y: number, camera: Camera | StereoCamera): AsyncPickData | undefined {
        this.prepare();

        if (this.dirty) {
            if (isTimingMode) this.webgl.timer.mark('PickHelper.asyncIdentify');
            this.render(camera);
            this.buffers.asyncRead();
            if (isTimingMode) this.webgl.timer.markEnd('PickHelper.asyncIdentify');
        }

        return {
            tryGet: () => {
                const status = this.buffers.check();
                if (status === AsyncPickStatus.Resolved) {
                    return this.getPickData(x, y, camera);
                } else if (status === AsyncPickStatus.Pending) {
                    return 'pending';
                } else if (status === AsyncPickStatus.Failed) {
                    this.dirty = true;
                }
            }
        };
    }

    reset() {
        this.buffers.reset();
        this.dirty = true;
    }

    dispose() {
        this.buffers.dispose();
    }

    constructor(private webgl: WebGLContext, private renderer: Renderer, private scene: Scene, private helper: Helper, private pickPass: PickPass, viewport: Viewport, options: PickOptions) {
        this.setViewport(viewport.x, viewport.y, viewport.width, viewport.height);
        this.pickPadding = options.pickPadding;

        if (!checkAsyncPickingSupport(webgl)) {
            this.asyncIdentify = (x, y, camera) => ({
                tryGet: () => this.identify(x, y, camera)
            });
        }
    }
}