/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BehaviorSubject, Subject, Subscription } from 'rxjs';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { Quat } from '../../mol-math/linear-algebra/3d/quat';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { Camera, ICamera } from '../camera';
import { PointerHelper } from './pointer-helper';
import { Vec2 } from '../../mol-math/linear-algebra/3d/vec2';
import { ButtonsType, InputObserver, ScreenTouchInput, TrackedPointerInput } from '../../mol-util/input/input-observer';
import { Plane3D } from '../../mol-math/geometry/primitives/plane3d';
import { Vec4 } from '../../mol-math/linear-algebra/3d/vec4';
import { StereoCamera } from '../camera/stereo';
import { Ray3D } from '../../mol-math/geometry/primitives/ray3d';
import { Scene } from '../../mol-gl/scene';
import { Sphere3D } from '../../mol-math/geometry';
import { Canvas3dInteractionHelper } from './interaction-events';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { cameraProject } from '../camera/util';
import { Binding } from '../../mol-util/binding';

const B = ButtonsType;
const Trigger = Binding.Trigger;
const Key = Binding.TriggerKey;

function getRigidTransformFromMat4(m: Mat4): XRRigidTransform {
    const d = Mat4.getDecomposition(m);
    return new XRRigidTransform(Vec3.toObj(d.position), Quat.toObj(d.quaternion));
}

function getRayFromPose(pose: XRPose, view?: Mat4): Ray3D {
    const origin = Vec3.fromObj(pose.transform.position);
    const t = Mat4.fromArray(Mat4(), pose.transform.matrix, 0);
    const td = Mat4.getDecomposition(t);
    const m = Mat4.fromQuat(Mat4(), td.quaternion);
    const direction = Vec3.transformMat4(Vec3(), Vec3.negUnitZ, m);
    const ray = Ray3D.create(origin, direction);
    if (view) Ray3D.transform(ray, ray, Mat4.invert(Mat4(), view));
    return ray;
}

type InputInfo = {
    targetRayPose: XRPose,
}

export const DefaultXRManagerBindings = {
    exit: Binding([Key('GamepadB')]),
    togglePassthrough: Binding([Key('GamepadA')]),
    gestureScale: Binding([Trigger(B.Flag.Trigger)]),
};
export const DefaultXRManagerAttribs = {
    bindings: DefaultXRManagerBindings,
};
export type XRManagerAttribs = typeof DefaultXRManagerAttribs

export const XRManagerParams = {
    minTargetDistance: PD.Numeric(0.4, { min: 0.001, max: 1, step: 0.001 }),
    disablePostprocessing: PD.Boolean(true),
    resolutionScale: PD.Numeric(1, { min: 0.1, max: 2, step: 0.1 }),
    sceneRadiusInMeters: PD.Numeric(0.25, { min: 0.01, max: 2, step: 0.01 }, { description: 'The radius of the scene bounding sphere in meters, used to set the initial camera scale.' }),
};
export type XRManagerParams = typeof XRManagerParams
export type XRManagerProps = PD.Values<XRManagerParams>

export class XRManager {
    private hoverSub: Subscription;
    private keyUpSub: Subscription;
    private gestureSub: Subscription;
    private sessionChangedSub: Subscription;

    readonly togglePassthrough = new Subject<void>();
    readonly sessionChanged = new Subject<void>();
    readonly isSupported = new BehaviorSubject(false);

    private xrSession: XRSession | undefined = undefined;
    get session() {
        return this.xrSession;
    }

    private xrRefSpace: XRReferenceSpace | undefined = undefined;

    private scaleFactor = 1;
    private prevScale = 0;
    private prevInput: { left?: InputInfo, right?: InputInfo } = {};
    private hit: Vec3 | undefined = undefined;

    readonly props: XRManagerProps;
    readonly attribs: XRManagerAttribs;

    setProps(props: Partial<XRManagerProps>) {
        Object.assign(this.props, props);
    }

    setAttribs(attribs: Partial<XRManagerAttribs>) {
        Object.assign(this.attribs, attribs);
    }

    private intersect(camera: ICamera, view: Mat4, plane: Plane3D, targetRayPose: XRPose): { point: Vec3, screen: Vec2 } | undefined {
        const point = Vec3();
        const ray = getRayFromPose(targetRayPose, view);
        if (Plane3D.intersectRay3D(point, plane, ray)) {
            const { height } = camera.viewport;
            const v = cameraProject(Vec4(), point, camera.viewport, camera.projectionView);
            const screen = Vec2.create(Math.floor(v[0]), height - Math.floor(v[1]));
            return { point, screen };
        }
    }

    setScaleFactor(factor: number) {
        this.scaleFactor = factor;
    }

    resetScale() {
        this.scaleFactor = 1;
        this.prevScale = 0;
    }

    update(xrFrame?: XRFrame): boolean {
        const { xrSession, xrRefSpace, input, camera, stereoCamera, pointerHelper } = this;
        if (!xrFrame || !xrSession || !xrRefSpace) return false;

        camera.scale = camera.scale * this.scaleFactor;
        this.prevScale = camera.scale;
        const camDirUnscaled = Vec3.sub(Vec3(), camera.position, camera.target);
        Vec3.scaleAndAdd(camera.position, camera.position, camDirUnscaled, 1 - this.scaleFactor);
        this.scaleFactor = 1;

        const xform = getRigidTransformFromMat4(camera.view);
        const xrOffsetRefSpace = xrRefSpace.getOffsetReferenceSpace(xform);
        const xrPose = xrFrame.getViewerPose(xrOffsetRefSpace);
        if (!xrPose) return false;

        const xrHeadPose = xrFrame.getViewerPose(xrRefSpace);
        if (xrHeadPose) {
            const hq = Quat.fromObj(xrHeadPose.transform.orientation);
            Mat4.fromQuat(camera.headRotation, hq);
        }

        const { depthFar, depthNear, baseLayer } = xrSession.renderState;
        if (!baseLayer) return false;

        if (depthFar !== camera.far || depthNear !== camera.near) {
            xrSession.updateRenderState({
                depthNear: camera.near,
                depthFar: camera.far,
            });
        }

        stereoCamera.update({ pose: xrPose, layer: baseLayer });
        const camLeft = stereoCamera.left;

        const cameraTarget = Vec3.scale(Vec3(), camLeft.state.target, camLeft.scale);
        const cameraPosition = Mat4.getTranslation(Vec3(), Mat4.invert(Mat4(), camLeft.view));
        const cameraDirection = Vec3.sub(Vec3(), cameraPosition, cameraTarget);
        const cameraPlane = Plane3D.fromNormalAndCoplanarPoint(Plane3D(), cameraDirection, cameraTarget);

        //

        const pointers: Ray3D[] = [];
        const points: Vec3[] = [];

        const trackedPointers: TrackedPointerInput[] = [];
        const screenTouches: ScreenTouchInput[] = [];

        if (xrSession.inputSources) {
            for (const inputSource of xrSession.inputSources) {
                if (inputSource.targetRayMode === 'screen') {
                    if (inputSource.gamepad) {
                        const { axes } = inputSource.gamepad;
                        const { width, height } = camLeft.viewport;
                        const x = ((axes[0] + 1) / 2) * width;
                        const y = ((axes[1] + 1) / 2) * height;
                        const ray = camLeft.getRay(Ray3D(), x, height - y);
                        screenTouches.push({ x, y, ray });
                    }
                    continue;
                }

                if (inputSource.targetRayMode !== 'tracked-pointer') continue;

                const { handedness, targetRaySpace, gamepad } = inputSource;
                if (!handedness) continue;

                const targetRayPose = xrFrame.getPose(targetRaySpace!, xrRefSpace);
                if (!targetRayPose) continue;

                const ray = getRayFromPose(targetRayPose, camera.view);
                pointers.push(ray);

                const sceneBoundingSphere = Sphere3D.scaleNX(Sphere3D(), this.scene.boundingSphereVisible, camLeft.scale);

                const si = Vec3();
                if (Ray3D.intersectSphere3D(si, ray, sceneBoundingSphere)) {
                    points.push(si);
                }

                let buttons = ButtonsType.create(ButtonsType.Flag.None);
                if (gamepad?.buttons[0]?.pressed) buttons |= ButtonsType.Flag.Primary;
                if (gamepad?.buttons[1]?.pressed) buttons |= ButtonsType.Flag.Secondary;
                if (gamepad?.buttons[3]?.pressed) buttons |= ButtonsType.Flag.Auxilary;
                if (gamepad?.buttons[4]?.pressed) buttons |= ButtonsType.Flag.Forth;
                if (gamepad?.buttons[5]?.pressed) buttons |= ButtonsType.Flag.Five;

                const prevInput = handedness === 'left' ? this.prevInput.left : this.prevInput.right;

                const intersection = this.intersect(camLeft, camera.view, cameraPlane, targetRayPose);
                const prevIntersection = prevInput ? this.intersect(camLeft, camera.view, cameraPlane, prevInput.targetRayPose) : undefined;

                const [x, y] = intersection?.screen ?? [0, 0];
                const [prevX, prevY] = prevIntersection?.screen ?? [x, y];

                const dd = Vec2.set(Vec2(), x - prevX, y - prevY);
                Vec2.setMagnitude(dd, dd, Math.min(100, Vec2.magnitude(dd)));
                const [dx, dy] = Vec2.round(dd, dd);

                trackedPointers.push({
                    handedness,
                    buttons,
                    x, y, dx, dy, ray,
                    axes: gamepad?.axes
                });

                if (handedness === 'left') {
                    this.prevInput.left = { targetRayPose };
                } else {
                    this.prevInput.right = { targetRayPose };
                }
            }
        } else {
            this.prevInput.left = undefined;
            this.prevInput.right = undefined;
        }

        input.updateTrackedPointers(trackedPointers);
        input.updateScreenTouches(screenTouches);

        pointerHelper.ensureEnabled();
        pointerHelper.update(pointers, points, this.hit);

        return true;
    }

    private async setSession(xrSession: XRSession | undefined) {
        if (this.xrSession === xrSession) return;

        await this.webgl.xr.set(xrSession, { resolutionScale: this.props.resolutionScale });

        this.xrSession = this.webgl.xr.session;
        this.prevInput = {};
        this.hit = undefined;

        if (this.xrSession) {
            this.xrRefSpace = await this.xrSession.requestReferenceSpace('local');
            this.pointerHelper.setProps({ enabled: 'on' });
            let scale = this.prevScale;
            if (scale === 0) {
                const { radius } = this.scene.boundingSphereVisible;
                scale = radius ? (1 / radius) * this.props.sceneRadiusInMeters : 0.01;
            }
            this.camera.forceFull = true;
            this.camera.scale = scale;
            this.camera.minTargetDistance = this.props.minTargetDistance;
            this.prevScale = scale;
        } else {
            this.xrRefSpace = undefined;
            Mat4.setZero(this.camera.headRotation);
            this.pointerHelper.setProps({ enabled: 'off' });
            this.camera.forceFull = false;
            this.camera.scale = 1;
            this.camera.minTargetDistance = 0;
        }
    }

    async end() {
        await this.webgl.xr.end();
    }

    private checkSupported = async () => {
        if (!navigator.xr) return false;

        const [arSupported, vrSupported] = await Promise.all([
            navigator.xr.isSessionSupported('immersive-ar'),
            navigator.xr.isSessionSupported('immersive-vr'),
        ]);
        this.isSupported.next(arSupported || vrSupported);
    };

    async request() {
        if (!navigator.xr) return;

        const session = await navigator.xr.isSessionSupported('immersive-ar')
            ? await navigator.xr.requestSession('immersive-ar')
            : await navigator.xr.requestSession('immersive-vr');

        await this.setSession(session);
    }

    dispose() {
        this.hoverSub.unsubscribe();
        this.keyUpSub.unsubscribe();
        this.gestureSub.unsubscribe();
        this.sessionChangedSub.unsubscribe();

        this.togglePassthrough.complete();
        this.sessionChanged.complete();
        this.isSupported.complete();

        navigator.xr?.removeEventListener('devicechange', this.checkSupported);
    }

    constructor(private webgl: WebGLContext, private input: InputObserver, private scene: Scene, private camera: Camera, private stereoCamera: StereoCamera, private pointerHelper: PointerHelper, private interactionHelper: Canvas3dInteractionHelper, props: Partial<XRManagerProps> = {}, attribs: Partial<XRManagerAttribs> = {}) {
        this.props = { ...PD.getDefaultValues(XRManagerParams), ...props };
        this.attribs = { ...DefaultXRManagerAttribs, ...attribs };

        this.hoverSub = this.interactionHelper.events.hover.subscribe(({ position }) => {
            this.hit = position;
        });

        this.sessionChangedSub = webgl.xr.changed.subscribe(async () => {
            await this.setSession(webgl.xr.session);
            this.sessionChanged.next();
        });

        this.checkSupported();
        navigator.xr?.addEventListener('devicechange', this.checkSupported);

        this.keyUpSub = input.keyUp.subscribe(({ code, modifiers, key }) => {
            const b = this.attribs.bindings;

            if (Binding.matchKey(b.exit, code, modifiers, key)) {
                this.end();
            }

            if (Binding.matchKey(b.togglePassthrough, code, modifiers, key)) {
                this.togglePassthrough.next();
            }
        });

        this.gestureSub = input.gesture.subscribe(({ scale, button, modifiers }) => {
            const b = this.attribs.bindings;

            if (Binding.match(b.gestureScale, button, modifiers)) {
                this.setScaleFactor(scale);
            }
        });
    }
}
