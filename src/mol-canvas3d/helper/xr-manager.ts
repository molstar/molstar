/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Subject, Subscription } from 'rxjs';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { Quat } from '../../mol-math/linear-algebra/3d/quat';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { Camera, ICamera } from '../camera';
import { PointerHelper } from './pointer-helper';
import { Vec2 } from '../../mol-math/linear-algebra/3d/vec2';
import { InputObserver } from '../../mol-util/input/input-observer';
import { Plane3D } from '../../mol-math/geometry/primitives/plane3d';
import { Vec4 } from '../../mol-math/linear-algebra/3d/vec4';
import { StereoCamera } from '../camera/stereo';
import { Ray3D } from '../../mol-math/geometry/primitives/ray3d';
import { Scene } from '../../mol-gl/scene';
import { Sphere3D } from '../../mol-math/geometry';
import { Canvas3dInteractionHelper } from './interaction-events';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { cameraProject } from '../camera/util';

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

export function getXRButton(manager: XRManager): HTMLElement | undefined {
    if (typeof document === 'undefined') return;
    if (!navigator.xr) return;

    const button = document.createElement('button');
    Object.assign(button.style, {
        bottom: '30px',
        left: '50%',
        transform: 'translateX(-50%)',
        background: 'black',
        zIndex: 1000,
        position: 'absolute',
        mixBlendMode: 'difference',
        color: 'white',
        border: 'white 2px solid',
        padding: '6px',
        display: 'none',
    });
    button.textContent = 'Enter XR';

    Promise.all([
        navigator.xr.isSessionSupported('immersive-ar'),
        navigator.xr.isSessionSupported('immersive-vr'),
    ]).then(([arSupported, vrSupported]) => {
        if (arSupported || vrSupported) {
            button.style.display = 'block';
        }
    });

    button.onclick = async () => {
        if (!navigator.xr) return;

        if (manager.session) {
            await manager.endSession();
        } else {
            const xrSession = await navigator.xr.isSessionSupported('immersive-ar')
                ? await navigator.xr.requestSession('immersive-ar')
                : await navigator.xr.requestSession('immersive-vr');
            console.log('WebXR session started', xrSession);
            xrSession.addEventListener('end', () => {
                button.textContent = 'Enter XR';
            });
            await manager.setSession(xrSession);
            if (manager.session) {
                button.textContent = 'Exit XR';
            }
        }
    };

    return button;
}

type InputInfo = {
    targetRayPose: XRPose,
    primaryButtonPressed: boolean,
    secondaryButtonPressed: boolean,
    aButtonPressed: boolean,
    bButtonPressed: boolean,
}

export const XRManagerParams = {
    enableHover: PD.Boolean(false),
    minTargetDistance: PD.Numeric(0.4, { min: 0.001, max: 1, step: 0.001 }),
    disablePostprocessing: PD.Boolean(true),
    resolutionScale: PD.Numeric(1, { min: 0.1, max: 2, step: 0.1 }),
};
export type XRManagerParams = typeof XRManagerParams
export type XRManagerProps = PD.Values<XRManagerParams>

export class XRManager {
    private hoverSub: Subscription;
    private sessionChangedSub: Subscription;

    readonly togglePassthrough = new Subject<void>();
    readonly sessionChanged = new Subject<void>();

    private xrSession: XRSession | undefined = undefined;
    get session() {
        return this.xrSession;
    }

    private xrRefSpace: XRReferenceSpace | undefined = undefined;

    private pointerDown = Vec2();
    private pointerEnd = Vec2();

    private prevInput: { left?: InputInfo, right?: InputInfo } = {};

    private hit: Vec3 | undefined = undefined;

    readonly props: XRManagerProps;

    setProps(props: Partial<XRManagerProps>) {
        Object.assign(this.props, props);
    }

    private intersect(camera: ICamera, view: Mat4, plane: Plane3D, input: InputInfo): { point: Vec3, screen: Vec2 } | undefined {
        const point = Vec3();
        const ray = getRayFromPose(input.targetRayPose, view);
        if (Plane3D.intersectRay3D(point, plane, ray)) {
            const { height } = camera.viewport;
            const v = cameraProject(Vec4(), point, camera.viewport, camera.projectionView);
            const screen = Vec2.create(Math.floor(v[0]), height - Math.floor(v[1]));
            return { point, screen };
        }
    }

    scaleFactor = 1;

    update(xrFrame?: XRFrame): boolean {
        const { xrSession: xrSession, xrRefSpace, input, camera, stereoCamera, pointerHelper } = this;
        if (!xrFrame || !xrSession || !xrRefSpace) return false;

        camera.setState({ scale: camera.state.scale * this.scaleFactor });
        this.prevScale = camera.state.scale;
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

        xrSession.updateRenderState({
            depthNear: camera.near,
            depthFar: camera.far,
        });
        const xrLayer = xrSession.renderState.baseLayer;
        if (!xrLayer) return false;

        stereoCamera.update({ pose: xrPose, layer: xrLayer });
        const camLeft = stereoCamera.left;

        const cameraTarget = Vec3.scale(Vec3(), camLeft.state.target, camLeft.state.scale);
        const cameraPosition = Mat4.getTranslation(Vec3(), Mat4.invert(Mat4(), camLeft.view));
        const cameraDirection = Vec3.sub(Vec3(), cameraPosition, cameraTarget);
        const cameraPlane = Plane3D.fromNormalAndCoplanarPoint(Plane3D(), cameraDirection, cameraTarget);

        //

        const pointers: Ray3D[] = [];
        const points: Vec3[] = [];
        for (const inputSource of xrSession.inputSources) {
            if (inputSource.targetRayMode !== 'tracked-pointer') continue;

            const targetRayPose = xrFrame.getPose(inputSource.targetRaySpace!, xrRefSpace);
            if (!targetRayPose) continue;

            const ray = getRayFromPose(targetRayPose, camera.view);
            pointers.push(ray);
            if (pointers.length > 1) continue;

            const sceneBoundingSphere = Sphere3D.scaleNX(Sphere3D(), this.scene.boundingSphereVisible, camLeft.state.scale);

            const si = Vec3();
            if (Ray3D.intersectSphere3D(si, ray, sceneBoundingSphere)) {
                points.push(si);
            }

            // if (inputSource.gamepad?.buttons[2].pressed) console.log(2);
            // if (inputSource.gamepad?.buttons[3].pressed) console.log(3);

            const primaryButton = inputSource.gamepad?.buttons[0];
            const secondaryButton = inputSource.gamepad?.buttons[1];
            const axesButton = inputSource.gamepad?.buttons[3];
            const aButton = inputSource.gamepad?.buttons[4];
            const bButton = inputSource.gamepad?.buttons[5];

            const inputInfo: InputInfo = {
                targetRayPose,
                primaryButtonPressed: !!primaryButton?.pressed,
                secondaryButtonPressed: !!secondaryButton?.pressed,
                aButtonPressed: !!aButton?.pressed,
                bButtonPressed: !!bButton?.pressed,
            };

            const intersection = this.intersect(camLeft, camera.view, cameraPlane, inputInfo);
            const prevIntersection = this.prevInput.right ? this.intersect(camLeft, camera.view, cameraPlane, this.prevInput.right) : undefined;

            if (primaryButton && intersection) {
                const [x, y] = intersection.screen;
                const [pageX, pageY] = intersection.screen;
                const modifiers = { alt: false, control: false, meta: false, shift: false };
                const button = primaryButton.pressed ? 1 : 0;
                const buttons = primaryButton.pressed ? 1 : 0;

                if (!!this.prevInput.right?.primaryButtonPressed && !primaryButton.pressed) {
                    Vec2.set(this.pointerEnd, x, y);
                    if (Vec2.distance(this.pointerEnd, this.pointerDown) < 10) {
                        input.click.next({ x, y, pageX, pageY, buttons: 1, button: 1, modifiers, ray });
                    }
                    input.interactionEnd.next(undefined);
                }

                if (this.props.enableHover || secondaryButton?.pressed) {
                    input.move.next({ x, y, pageX, pageY, buttons, button, modifiers, inside: true, onElement: true, ray });
                } else if (!secondaryButton?.pressed && this.prevInput.right?.secondaryButtonPressed) {
                    input.leave.next(undefined);
                }

                if (primaryButton.pressed) {
                    const isStart = !this.prevInput.right?.primaryButtonPressed;
                    const [prevX, prevY] = prevIntersection?.screen ?? [x, y];

                    const dd = Vec2.set(Vec2(), x - prevX, y - prevY);
                    // Vec2.scale(dd, dd, 1 / (camera.state.scale * 100));
                    const dm = Vec2.magnitude(dd);
                    Vec2.setMagnitude(dd, dd, Math.min(100, dm));
                    Vec2.round(dd, dd);
                    const [dx, dy] = dd;

                    input.drag.next({ x, y, dx, dy, pageX, pageY, buttons, button, modifiers, isStart, useDelta: true });
                    // console.log('xr input', { x, y, dx, dy, pageX, pageY, buttons, button, modifiers });
                    if (isStart) {
                        Vec2.set(this.pointerDown, x, y);
                    }
                }
            }

            if (intersection) points.push(intersection.point);
            if (prevIntersection) points.push(prevIntersection.point);

            if (this.prevInput.right?.aButtonPressed && !aButton?.pressed) {
                this.togglePassthrough.next();
            }

            if (this.prevInput.right?.bButtonPressed && !bButton?.pressed) {
                this.endSession();
                return false;
            }

            if (inputSource.gamepad?.axes) {
                const x = 0, y = 0, pageX = 0, pageY = 0;
                const preventDefault = () => {};
                const modifiers = { alt: false, control: false, meta: false, shift: false };
                const keyValue = { x, y, pageX, pageY, modifiers, preventDefault };
                const { axes } = inputSource.gamepad;

                if (axesButton?.pressed) {
                    if (axes[3] < 0) {
                        this.scaleFactor = 1.01;
                    } else if (axes[3] > 0) {
                        this.scaleFactor = 0.99;
                    }
                } else {
                    if (axes[3] < 0) {
                        input.keyDown.next({ ...keyValue, code: 'KeyW', key: 'w' });
                        input.keyUp.next({ ...keyValue, code: 'KeyS', key: 's' });
                    } else if (axes[3] > 0) {
                        input.keyDown.next({ ...keyValue, code: 'KeyS', key: 's' });
                        input.keyUp.next({ ...keyValue, code: 'KeyW', key: 'w' });
                    } else {
                        input.keyUp.next({ ...keyValue, code: 'KeyW', key: 'w' });
                        input.keyUp.next({ ...keyValue, code: 'KeyS', key: 's' });
                    }
                }
            }

            this.prevInput.right = inputInfo;
        }

        const plane = { point: cameraTarget, normal: cameraPlane.normal };

        // console.log('pointers', pointers, points);
        pointerHelper.update(pointers, points, this.hit, plane);

        return true;
    }

    prevScale = 0;

    async setSession(xrSession: XRSession | undefined) {
        if (this.xrSession === xrSession) return;

        await this.webgl.setXRSession(xrSession, { resolutionScale: this.props.resolutionScale });

        this.xrSession = this.webgl.xrSession;
        this.prevInput = {};
        this.hit = undefined;

        if (this.xrSession) {
            this.xrRefSpace = await this.xrSession.requestReferenceSpace('local');
            // console.log('xrRefSpace', this.xrRefSpace);
            this.pointerHelper.setProps({ enabled: 'on' });
            let scale = this.prevScale;
            if (scale === 0) {
                const { radius } = this.scene.boundingSphereVisible;
                scale = radius ? 1 / (radius * 3) : 0.01;
            }
            this.camera.setState({
                forceFull: true,
                scale,
                minTargetDistance: this.props.minTargetDistance,
            });
            this.prevScale = scale;
        } else {
            this.xrRefSpace = undefined;
            Mat4.setZero(this.camera.headRotation);
            this.pointerHelper.setProps({ enabled: 'off' });
            this.camera.setState({
                forceFull: false,
                scale: 1,
                minTargetDistance: 0,
            });
        }
    }

    async endSession() {
        await this.webgl.endXRSession();
    }

    dispose() {
        this.hoverSub.unsubscribe();
        this.sessionChangedSub.unsubscribe();

        this.togglePassthrough.complete();
        this.sessionChanged.complete();
    }

    constructor(private webgl: WebGLContext, private input: InputObserver, private scene: Scene, private camera: Camera, private stereoCamera: StereoCamera, private pointerHelper: PointerHelper, private interactionHelper: Canvas3dInteractionHelper, props: Partial<XRManagerProps> = {}) {
        this.props = { ...PD.getDefaultValues(XRManagerParams), ...props };

        this.hoverSub = this.interactionHelper.events.hover.subscribe(({ position }) => {
            this.hit = position;
        });

        this.sessionChangedSub = webgl.xrSessionChanged.subscribe(async () => {
            await this.setSession(webgl.xrSession);
            this.sessionChanged.next();
        });
    }
}
