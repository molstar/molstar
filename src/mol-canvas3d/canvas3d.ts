/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, Subscription } from 'rxjs';
import { now } from '../mol-util/now';
import { Vec3, Vec2 } from '../mol-math/linear-algebra';
import { InputObserver, ModifiersKeys, ButtonsType } from '../mol-util/input/input-observer';
import { Renderer, RendererStats, RendererParams } from '../mol-gl/renderer';
import { GraphicsRenderObject } from '../mol-gl/render-object';
import { TrackballControls, TrackballControlsParams } from './controls/trackball';
import { Viewport } from './camera/util';
import { createContext, WebGLContext, getGLContext } from '../mol-gl/webgl/context';
import { Representation } from '../mol-repr/representation';
import { Scene } from '../mol-gl/scene';
import { PickingId } from '../mol-geo/geometry/picking';
import { MarkerAction } from '../mol-util/marker-action';
import { Loci, EmptyLoci, isEmptyLoci } from '../mol-model/loci';
import { Camera } from './camera';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { DebugHelperParams } from './helper/bounding-sphere-helper';
import { SetUtils } from '../mol-util/set';
import { Canvas3dInteractionHelper, Canvas3dInteractionHelperParams } from './helper/interaction-events';
import { PostprocessingParams } from './passes/postprocessing';
import { MultiSampleHelper, MultiSampleParams, MultiSamplePass } from './passes/multi-sample';
import { PickData } from './passes/pick';
import { PickHelper } from './passes/pick';
import { ImagePass, ImageProps } from './passes/image';
import { Sphere3D } from '../mol-math/geometry';
import { isDebugMode } from '../mol-util/debug';
import { CameraHelperParams } from './helper/camera-helper';
import { produce } from 'immer';
import { HandleHelperParams } from './helper/handle-helper';
import { StereoCamera, StereoCameraParams } from './camera/stereo';
import { Helper } from './helper/helper';
import { Passes } from './passes/passes';
import { shallowEqual } from '../mol-util';
import { MarkingParams } from './passes/marking';

export const Canvas3DParams = {
    camera: PD.Group({
        mode: PD.Select('perspective', PD.arrayToOptions(['perspective', 'orthographic'] as const), { label: 'Camera' }),
        helper: PD.Group(CameraHelperParams, { isFlat: true }),
        stereo: PD.MappedStatic('off', {
            on: PD.Group(StereoCameraParams),
            off: PD.Group({})
        }, { cycle: true, hideIf: p => p?.mode !== 'perspective' }),
        manualReset: PD.Boolean(false, { isHidden: true }),
    }, { pivot: 'mode' }),
    cameraFog: PD.MappedStatic('on', {
        on: PD.Group({
            intensity: PD.Numeric(15, { min: 1, max: 100, step: 1 }),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Show fog in the distance' }),
    cameraClipping: PD.Group({
        radius: PD.Numeric(100, { min: 0, max: 99, step: 1 }, { label: 'Clipping', description: 'How much of the scene to show.' }),
        far: PD.Boolean(true, { description: 'Hide scene in the distance' }),
    }, { pivot: 'radius' }),
    viewport: PD.MappedStatic('canvas', {
        canvas: PD.Group({}),
        'static-frame': PD.Group({
            x: PD.Numeric(0),
            y: PD.Numeric(0),
            width: PD.Numeric(128),
            height: PD.Numeric(128)
        }),
        'relative-frame': PD.Group({
            x: PD.Numeric(0.33, { min: 0, max: 1, step: 0.01 }),
            y: PD.Numeric(0.33, { min: 0, max: 1, step: 0.01 }),
            width: PD.Numeric(0.5, { min: 0.01, max: 1, step: 0.01 }),
            height: PD.Numeric(0.5, { min: 0.01, max: 1, step: 0.01 })
        })
    }),

    cameraResetDurationMs: PD.Numeric(250, { min: 0, max: 1000, step: 1 }, { description: 'The time it takes to reset the camera.' }),
    transparentBackground: PD.Boolean(false),

    multiSample: PD.Group(MultiSampleParams),
    postprocessing: PD.Group(PostprocessingParams),
    marking: PD.Group(MarkingParams),
    renderer: PD.Group(RendererParams),
    trackball: PD.Group(TrackballControlsParams),
    interaction: PD.Group(Canvas3dInteractionHelperParams),
    debug: PD.Group(DebugHelperParams),
    handle: PD.Group(HandleHelperParams),
};
export const DefaultCanvas3DParams = PD.getDefaultValues(Canvas3DParams);
export type Canvas3DProps = PD.Values<typeof Canvas3DParams>
export type PartialCanvas3DProps = {
    [K in keyof Canvas3DProps]?: Canvas3DProps[K] extends { name: string, params: any } ? Canvas3DProps[K] : Partial<Canvas3DProps[K]>
}

export { Canvas3DContext };

/** Can be used to create multiple Canvas3D objects */
interface Canvas3DContext {
    readonly canvas: HTMLCanvasElement
    readonly webgl: WebGLContext
    readonly input: InputObserver
    readonly passes: Passes
    readonly attribs: Readonly<Canvas3DContext.Attribs>
    readonly contextLost: BehaviorSubject<now.Timestamp>
    readonly contextRestored: BehaviorSubject<now.Timestamp>
    dispose: (options?: Partial<{ doNotForceWebGLContextLoss: boolean }>) => void
}

namespace Canvas3DContext {
    export const DefaultAttribs = {
        /** true by default to avoid issues with Safari (Jan 2021) */
        antialias: true,
        /** true to support multiple Canvas3D objects with a single context */
        preserveDrawingBuffer: true,
        pixelScale: 1,
        pickScale: 0.25,
        /** extra pixels to around target to check in case target is empty */
        pickPadding: 1,
        enableWboit: true,
        preferWebGl1: false
    };
    export type Attribs = typeof DefaultAttribs

    export function fromCanvas(canvas: HTMLCanvasElement, attribs: Partial<Attribs> = {}): Canvas3DContext {
        const a = { ...DefaultAttribs, ...attribs };
        const { antialias, preserveDrawingBuffer, pixelScale, preferWebGl1 } = a;
        const gl = getGLContext(canvas, {
            antialias,
            preserveDrawingBuffer,
            alpha: true, // the renderer requires an alpha channel
            depth: true, // the renderer requires a depth buffer
            premultipliedAlpha: true, // the renderer outputs PMA
            preferWebGl1
        });
        if (gl === null) throw new Error('Could not create a WebGL rendering context');

        const input = InputObserver.fromElement(canvas, { pixelScale, preventGestures: true });
        const webgl = createContext(gl, { pixelScale });
        const passes = new Passes(webgl, attribs);

        if (isDebugMode) {
            const loseContextExt = gl.getExtension('WEBGL_lose_context');
            if (loseContextExt) {
                // Hold down shift+ctrl+alt and press any mouse button to call `loseContext`.
                // After 1 second `restoreContext` will be called.
                canvas.addEventListener('mousedown', e => {
                    if (webgl.isContextLost) return;
                    if (!e.shiftKey || !e.ctrlKey || !e.altKey) return;

                    if (isDebugMode) console.log('lose context');
                    loseContextExt.loseContext();

                    setTimeout(() => {
                        if (!webgl.isContextLost) return;
                        if (isDebugMode) console.log('restore context');
                        loseContextExt.restoreContext();
                    }, 1000);
                }, false);
            }
        }

        // https://www.khronos.org/webgl/wiki/HandlingContextLost

        const contextLost = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp);

        const handleWebglContextLost = (e: Event) => {
            webgl.setContextLost();
            e.preventDefault();
            if (isDebugMode) console.log('context lost');
            contextLost.next(now());
        };

        const handlewWebglContextRestored = () => {
            if (!webgl.isContextLost) return;
            webgl.handleContextRestored(() => {
                passes.draw.reset();
            });
            if (isDebugMode) console.log('context restored');
        };

        canvas.addEventListener('webglcontextlost', handleWebglContextLost, false);
        canvas.addEventListener('webglcontextrestored', handlewWebglContextRestored, false);

        return {
            canvas,
            webgl,
            input,
            passes,
            attribs: a,
            contextLost,
            contextRestored: webgl.contextRestored,
            dispose: (options?: Partial<{ doNotForceWebGLContextLoss: boolean }>) => {
                input.dispose();

                canvas.removeEventListener('webglcontextlost', handleWebglContextLost, false);
                canvas.removeEventListener('webglcontextrestored', handlewWebglContextRestored, false);
                webgl.destroy(options);
            }
        };
    }
}

export { Canvas3D };

interface Canvas3D {
    readonly webgl: WebGLContext,

    add(repr: Representation.Any): void
    remove(repr: Representation.Any): void
    /**
     * This function must be called if animate() is not set up so that add/remove actions take place.
     */
    commit(isSynchronous?: boolean): void
    /**
     * Function for external "animation" control
     * Calls commit.
     */
    tick(t: now.Timestamp, options?: { isSynchronous?: boolean, manualDraw?: boolean }): void
    update(repr?: Representation.Any, keepBoundingSphere?: boolean): void
    clear(): void
    syncVisibility(): void

    requestDraw(force?: boolean): void

    /** Reset the timers, used by "animate" */
    resetTime(t: number): void
    animate(): void
    /**
     * Pause animation loop and optionally any rendering
     * @param noDraw pause any rendering (drawPaused = true)
     */
    pause(noDraw?: boolean): void
    /** Sets drawPaused = false without starting the built in animation loop */
    resume(): void
    identify(x: number, y: number): PickData | undefined
    mark(loci: Representation.Loci, action: MarkerAction, noDraw?: boolean): void
    getLoci(pickingId: PickingId | undefined): Representation.Loci

    notifyDidDraw: boolean,
    readonly didDraw: BehaviorSubject<now.Timestamp>
    readonly commited: BehaviorSubject<now.Timestamp>
    readonly reprCount: BehaviorSubject<number>
    readonly resized: BehaviorSubject<any>

    handleResize(): void
    /** performs handleResize on the next animation frame */
    requestResize(): void
    /** Focuses camera on scene's bounding sphere, centered and zoomed. */
    requestCameraReset(options?: { durationMs?: number, snapshot?: Camera.SnapshotProvider }): void
    readonly camera: Camera
    readonly boundingSphere: Readonly<Sphere3D>
    readonly boundingSphereVisible: Readonly<Sphere3D>
    setProps(props: PartialCanvas3DProps | ((old: Canvas3DProps) => Partial<Canvas3DProps> | void), doNotRequestDraw?: boolean /* = false */): void
    getImagePass(props: Partial<ImageProps>): ImagePass
    getRenderObjects(): GraphicsRenderObject[]

    /** Returns a copy of the current Canvas3D instance props */
    readonly props: Readonly<Canvas3DProps>
    readonly input: InputObserver
    readonly stats: RendererStats
    readonly interaction: Canvas3dInteractionHelper['events']

    dispose(): void
}

const requestAnimationFrame = typeof window !== 'undefined'
    ? window.requestAnimationFrame
    : (f: (time: number) => void) => setImmediate(() => f(Date.now())) as unknown as number;
const cancelAnimationFrame = typeof window !== 'undefined'
    ? window.cancelAnimationFrame
    : (handle: number) => clearImmediate(handle as unknown as NodeJS.Immediate);

namespace Canvas3D {
    export interface HoverEvent { current: Representation.Loci, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys, page?: Vec2, position?: Vec3 }
    export interface DragEvent { current: Representation.Loci, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys, pageStart: Vec2, pageEnd: Vec2 }
    export interface ClickEvent { current: Representation.Loci, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys, page?: Vec2, position?: Vec3 }

    export function create({ webgl, input, passes, attribs }: Canvas3DContext, props: Partial<Canvas3DProps> = {}): Canvas3D {
        const p: Canvas3DProps = { ...DefaultCanvas3DParams, ...props };

        const reprRenderObjects = new Map<Representation.Any, Set<GraphicsRenderObject>>();
        const reprUpdatedSubscriptions = new Map<Representation.Any, Subscription>();
        const reprCount = new BehaviorSubject(0);

        let startTime = now();
        const didDraw = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp);
        const commited = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp);

        const { gl, contextRestored } = webgl;

        let x = 0;
        let y = 0;
        let width = 128;
        let height = 128;
        updateViewport();

        const scene = Scene.create(webgl);

        const camera = new Camera({
            position: Vec3.create(0, 0, 100),
            mode: p.camera.mode,
            fog: p.cameraFog.name === 'on' ? p.cameraFog.params.intensity : 0,
            clipFar: p.cameraClipping.far
        }, { x, y, width, height }, { pixelScale: attribs.pixelScale });
        const stereoCamera = new StereoCamera(camera, p.camera.stereo.params);

        const controls = TrackballControls.create(input, camera, p.trackball);
        const renderer = Renderer.create(webgl, p.renderer);
        const helper = new Helper(webgl, scene, p);

        const pickHelper = new PickHelper(webgl, renderer, scene, helper, passes.pick, { x, y, width, height }, attribs.pickPadding);
        const interactionHelper = new Canvas3dInteractionHelper(identify, getLoci, input, camera, p.interaction);
        const multiSampleHelper = new MultiSampleHelper(passes.multiSample);

        let cameraResetRequested = false;
        let nextCameraResetDuration: number | undefined = void 0;
        let nextCameraResetSnapshot: Camera.SnapshotProvider | undefined = void 0;
        let resizeRequested = false;

        let notifyDidDraw = true;

        function getLoci(pickingId: PickingId | undefined) {
            let loci: Loci = EmptyLoci;
            let repr: Representation.Any = Representation.Empty;
            if (pickingId) {
                const cameraHelperLoci = helper.camera.getLoci(pickingId);
                if (cameraHelperLoci !== EmptyLoci) return { loci: cameraHelperLoci, repr };

                loci = helper.handle.getLoci(pickingId);

                reprRenderObjects.forEach((_, _repr) => {
                    const _loci = _repr.getLoci(pickingId);
                    if (!isEmptyLoci(_loci)) {
                        if (!isEmptyLoci(loci)) {
                            console.warn('found another loci, this should not happen');
                        }
                        loci = _loci;
                        repr = _repr;
                    }
                });
            }
            return { loci, repr };
        }

        function mark(reprLoci: Representation.Loci, action: MarkerAction, noDraw = false) {
            const { repr, loci } = reprLoci;
            let changed = false;
            if (repr) {
                changed = repr.mark(loci, action);
            } else {
                changed = helper.handle.mark(loci, action);
                changed = helper.camera.mark(loci, action) || changed;
                reprRenderObjects.forEach((_, _repr) => { changed = _repr.mark(loci, action) || changed; });
            }
            if (changed && !noDraw) {
                scene.update(void 0, true);
                helper.handle.scene.update(void 0, true);
                helper.camera.scene.update(void 0, true);
                const prevPickDirty = pickHelper.dirty;
                draw(true);
                pickHelper.dirty = prevPickDirty; // marking does not change picking buffers
            }
        }

        function render(force: boolean) {
            if (webgl.isContextLost) return false;

            let resized = false;
            if (resizeRequested) {
                handleResize(false);
                resizeRequested = false;
                resized = true;
            }

            if (x > gl.drawingBufferWidth || x + width < 0 ||
                y > gl.drawingBufferHeight || y + height < 0
            ) return false;

            let didRender = false;
            controls.update(currentTime);
            const cameraChanged = camera.update();

            const shouldRender = force || cameraChanged || resized || forceNextRender;
            forceNextRender = false;

            const multiSampleChanged = multiSampleHelper.update(shouldRender, p.multiSample);

            if (shouldRender || multiSampleChanged) {
                let cam: Camera | StereoCamera = camera;
                if (p.camera.stereo.name === 'on') {
                    stereoCamera.update();
                    cam = stereoCamera;
                }

                if (MultiSamplePass.isEnabled(p.multiSample)) {
                    multiSampleHelper.render(renderer, cam, scene, helper, true, p.transparentBackground, p);
                } else {
                    passes.draw.render(renderer, cam, scene, helper, true, p.transparentBackground, p.postprocessing, p.marking);
                }
                pickHelper.dirty = true;
                didRender = true;
            }

            return didRender;
        }

        let forceNextRender = false;
        let forceDrawAfterAllCommited = false;
        let currentTime = 0;
        let drawPaused = false;

        function draw(force?: boolean) {
            if (drawPaused) return;
            if (render(!!force) && notifyDidDraw) {
                didDraw.next(now() - startTime as now.Timestamp);
            }
        }

        function requestDraw(force?: boolean) {
            forceNextRender = forceNextRender || !!force;
        }

        let animationFrameHandle = 0;

        function tick(t: now.Timestamp, options?: { isSynchronous?: boolean, manualDraw?: boolean }) {
            currentTime = t;
            commit(options?.isSynchronous);
            camera.transition.tick(currentTime);

            if (options?.manualDraw) {
                return;
            }

            draw(false);
            if (!camera.transition.inTransition && !webgl.isContextLost) {
                interactionHelper.tick(currentTime);
            }
        }

        function _animate() {
            tick(now());
            animationFrameHandle = requestAnimationFrame(_animate);
        }

        function resetTime(t: now.Timestamp) {
            startTime = t;
            controls.start(t);
        }

        function animate() {
            drawPaused = false;
            controls.start(now());
            if (animationFrameHandle === 0) _animate();
        }

        function pause(noDraw = false) {
            drawPaused = noDraw;
            cancelAnimationFrame(animationFrameHandle);
            animationFrameHandle = 0;
        }

        function identify(x: number, y: number): PickData | undefined {
            const cam = p.camera.stereo.name === 'on' ? stereoCamera : camera;
            return webgl.isContextLost ? undefined : pickHelper.identify(x, y, cam);
        }

        function commit(isSynchronous: boolean = false) {
            const allCommited = commitScene(isSynchronous);
            // Only reset the camera after the full scene has been commited.
            if (allCommited) {
                resolveCameraReset();
                if (forceDrawAfterAllCommited) {
                    if (helper.debug.isEnabled) helper.debug.update();
                    draw(true);
                    forceDrawAfterAllCommited = false;
                }
                commited.next(now());
            }
        }

        function resolveCameraReset() {
            if (!cameraResetRequested) return;

            const boundingSphere = scene.boundingSphereVisible;
            const { center, radius } = boundingSphere;

            const autoAdjustControls = controls.props.autoAdjustMinMaxDistance;
            if (autoAdjustControls.name === 'on') {
                const minDistance = autoAdjustControls.params.minDistanceFactor * radius + autoAdjustControls.params.minDistancePadding;
                const maxDistance = Math.max(autoAdjustControls.params.maxDistanceFactor * radius, autoAdjustControls.params.maxDistanceMin);
                controls.setProps({ minDistance, maxDistance });
            }

            if (radius > 0) {
                const duration = nextCameraResetDuration === undefined ? p.cameraResetDurationMs : nextCameraResetDuration;
                const focus = camera.getFocus(center, radius);
                const next = typeof nextCameraResetSnapshot === 'function' ? nextCameraResetSnapshot(scene, camera) : nextCameraResetSnapshot;
                const snapshot = next ? { ...focus, ...next } : focus;
                camera.setState({ ...snapshot, radiusMax: scene.boundingSphere.radius }, duration);
            }

            nextCameraResetDuration = void 0;
            nextCameraResetSnapshot = void 0;
            cameraResetRequested = false;
        }

        const oldBoundingSphereVisible = Sphere3D();
        const cameraSphere = Sphere3D();

        function shouldResetCamera() {
            if (camera.state.radiusMax === 0) return true;

            if (camera.transition.inTransition || nextCameraResetSnapshot) return false;

            let cameraSphereOverlapsNone = true;
            Sphere3D.set(cameraSphere, camera.state.target, camera.state.radius);

            // check if any renderable has moved outside of the old bounding sphere
            // and if no renderable is overlapping with the camera sphere
            for (const r of scene.renderables) {
                if (!r.state.visible) continue;

                const b = r.values.boundingSphere.ref.value;
                if (!b.radius) continue;

                const cameraDist = Vec3.distance(cameraSphere.center, b.center);
                if ((cameraDist > cameraSphere.radius || cameraDist > b.radius || b.radius > camera.state.radiusMax) && !Sphere3D.includes(oldBoundingSphereVisible, b)) return true;
                if (Sphere3D.overlaps(cameraSphere, b)) cameraSphereOverlapsNone = false;
            }

            return cameraSphereOverlapsNone;
        }

        const sceneCommitTimeoutMs = 250;
        function commitScene(isSynchronous: boolean) {
            if (!scene.needsCommit) return true;

            // snapshot the current bounding sphere of visible objects
            Sphere3D.copy(oldBoundingSphereVisible, scene.boundingSphereVisible);

            if (!scene.commit(isSynchronous ? void 0 : sceneCommitTimeoutMs)) return false;

            if (helper.debug.isEnabled) helper.debug.update();
            if (!p.camera.manualReset && (reprCount.value === 0 || shouldResetCamera())) {
                cameraResetRequested = true;
            }
            if (oldBoundingSphereVisible.radius === 0) nextCameraResetDuration = 0;

            camera.setState({ radiusMax: scene.boundingSphere.radius }, 0);
            reprCount.next(reprRenderObjects.size);
            if (isDebugMode) consoleStats();

            return true;
        }

        function consoleStats() {
            console.table(scene.renderables.map(r => ({
                drawCount: r.values.drawCount.ref.value,
                instanceCount: r.values.instanceCount.ref.value,
                materialId: r.materialId,
                renderItemId: r.id,
            })));
            console.log(webgl.stats);

            const { texture, attribute, elements } = webgl.resources.getByteCounts();
            console.log({
                texture: `${(texture / 1024 / 1024).toFixed(3)} MiB`,
                attribute: `${(attribute / 1024 / 1024).toFixed(3)} MiB`,
                elements: `${(elements / 1024 / 1024).toFixed(3)} MiB`,
            });
        }

        function add(repr: Representation.Any) {
            registerAutoUpdate(repr);

            const oldRO = reprRenderObjects.get(repr);
            const newRO = new Set<GraphicsRenderObject>();
            repr.renderObjects.forEach(o => newRO.add(o));

            if (oldRO) {
                if (!SetUtils.areEqual(newRO, oldRO)) {
                    newRO.forEach(o => { if (!oldRO.has(o)) scene.add(o); });
                    oldRO.forEach(o => { if (!newRO.has(o)) scene.remove(o); });
                }
            } else {
                repr.renderObjects.forEach(o => scene.add(o));
            }
            reprRenderObjects.set(repr, newRO);

            scene.update(repr.renderObjects, false);
            forceDrawAfterAllCommited = true;
            if (isDebugMode) consoleStats();
        }

        function remove(repr: Representation.Any) {
            unregisterAutoUpdate(repr);

            const renderObjects = reprRenderObjects.get(repr);
            if (renderObjects) {
                renderObjects.forEach(o => scene.remove(o));
                reprRenderObjects.delete(repr);
                forceDrawAfterAllCommited = true;
                if (isDebugMode) consoleStats();
            }
        }

        function registerAutoUpdate(repr: Representation.Any) {
            if (reprUpdatedSubscriptions.has(repr)) return;

            reprUpdatedSubscriptions.set(repr, repr.updated.subscribe(_ => {
                if (!repr.state.syncManually) add(repr);
            }));
        }

        function unregisterAutoUpdate(repr: Representation.Any) {
            const updatedSubscription = reprUpdatedSubscriptions.get(repr);
            if (updatedSubscription) {
                updatedSubscription.unsubscribe();
                reprUpdatedSubscriptions.delete(repr);
            }
        }

        function getProps(): Canvas3DProps {
            const radius = scene.boundingSphere.radius > 0
                ? 100 - Math.round((camera.transition.target.radius / scene.boundingSphere.radius) * 100)
                : 0;

            return {
                camera: {
                    mode: camera.state.mode,
                    helper: { ...helper.camera.props },
                    stereo: { ...p.camera.stereo },
                    manualReset: !!p.camera.manualReset
                },
                cameraFog: camera.state.fog > 0
                    ? { name: 'on' as const, params: { intensity: camera.state.fog } }
                    : { name: 'off' as const, params: {} },
                cameraClipping: { far: camera.state.clipFar, radius },
                cameraResetDurationMs: p.cameraResetDurationMs,
                transparentBackground: p.transparentBackground,
                viewport: p.viewport,

                postprocessing: { ...p.postprocessing },
                marking: { ...p.marking },
                multiSample: { ...p.multiSample },
                renderer: { ...renderer.props },
                trackball: { ...controls.props },
                interaction: { ...interactionHelper.props },
                debug: { ...helper.debug.props },
                handle: { ...helper.handle.props },
            };
        }

        const contextRestoredSub = contextRestored.subscribe(() => {
            pickHelper.dirty = true;
            draw(true);
            // Unclear why, but in Chrome with wboit enabled the first `draw` only clears
            // the drawingBuffer. Note that in Firefox the drawingBuffer is preserved after
            // context loss so it is unclear if it behaves the same.
            draw(true);
        });

        const resized = new BehaviorSubject<any>(0);

        function handleResize(draw = true) {
            passes.updateSize();
            updateViewport();
            syncViewport();
            if (draw) requestDraw(true);
            resized.next(+new Date());
        }

        return {
            webgl,

            add,
            remove,
            commit,
            update: (repr, keepSphere) => {
                if (repr) {
                    if (!reprRenderObjects.has(repr)) return;
                    scene.update(repr.renderObjects, !!keepSphere);
                } else {
                    scene.update(void 0, !!keepSphere);
                }
                forceDrawAfterAllCommited = true;
            },
            clear: () => {
                reprUpdatedSubscriptions.forEach(v => v.unsubscribe());
                reprUpdatedSubscriptions.clear();
                reprRenderObjects.clear();
                scene.clear();
                helper.debug.clear();
                requestDraw(true);
                reprCount.next(reprRenderObjects.size);
            },
            syncVisibility: () => {
                if (camera.state.radiusMax === 0) {
                    cameraResetRequested = true;
                    nextCameraResetDuration = 0;
                }

                if (scene.syncVisibility()) {
                    if (helper.debug.isEnabled) helper.debug.update();
                }
                requestDraw(true);
            },

            requestDraw,
            tick,
            animate,
            resetTime,
            pause,
            resume: () => { drawPaused = false; },
            identify,
            mark,
            getLoci,

            handleResize,
            requestResize: () => {
                resizeRequested = true;
            },
            requestCameraReset: options => {
                nextCameraResetDuration = options?.durationMs;
                nextCameraResetSnapshot = options?.snapshot;
                cameraResetRequested = true;
            },
            camera,
            boundingSphere: scene.boundingSphere,
            boundingSphereVisible: scene.boundingSphereVisible,
            get notifyDidDraw() { return notifyDidDraw; },
            set notifyDidDraw(v: boolean) { notifyDidDraw = v; },
            didDraw,
            commited,
            reprCount,
            resized,
            setProps: (properties, doNotRequestDraw = false) => {
                const props: PartialCanvas3DProps = typeof properties === 'function'
                    ? produce(getProps(), properties as any)
                    : properties;

                const cameraState: Partial<Camera.Snapshot> = Object.create(null);
                if (props.camera && props.camera.mode !== undefined && props.camera.mode !== camera.state.mode) {
                    cameraState.mode = props.camera.mode;
                }
                if (props.cameraFog !== undefined && props.cameraFog.params) {
                    const newFog = props.cameraFog.name === 'on' ? props.cameraFog.params.intensity : 0;
                    if (newFog !== camera.state.fog) cameraState.fog = newFog;
                }
                if (props.cameraClipping !== undefined) {
                    if (props.cameraClipping.far !== undefined && props.cameraClipping.far !== camera.state.clipFar) {
                        cameraState.clipFar = props.cameraClipping.far;
                    }
                    if (props.cameraClipping.radius !== undefined) {
                        const radius = (scene.boundingSphere.radius / 100) * (100 - props.cameraClipping.radius);
                        if (radius > 0 && radius !== cameraState.radius) {
                            // if radius = 0, NaNs happen
                            cameraState.radius = Math.max(radius, 0.01);
                        }
                    }
                }
                if (Object.keys(cameraState).length > 0) camera.setState(cameraState);

                if (props.camera?.helper) helper.camera.setProps(props.camera.helper);
                if (props.camera?.manualReset !== undefined) p.camera.manualReset = props.camera.manualReset;
                if (props.camera?.stereo !== undefined) Object.assign(p.camera.stereo, props.camera.stereo);
                if (props.cameraResetDurationMs !== undefined) p.cameraResetDurationMs = props.cameraResetDurationMs;
                if (props.transparentBackground !== undefined) p.transparentBackground = props.transparentBackground;
                if (props.viewport !== undefined) {
                    const doNotUpdate = p.viewport === props.viewport ||
                        (p.viewport.name === props.viewport.name && shallowEqual(p.viewport.params, props.viewport.params));

                    if (!doNotUpdate) {
                        p.viewport = props.viewport;
                        updateViewport();
                        syncViewport();
                    }
                }

                if (props.postprocessing) Object.assign(p.postprocessing, props.postprocessing);
                if (props.marking) Object.assign(p.marking, props.marking);
                if (props.multiSample) Object.assign(p.multiSample, props.multiSample);
                if (props.renderer) renderer.setProps(props.renderer);
                if (props.trackball) controls.setProps(props.trackball);
                if (props.interaction) interactionHelper.setProps(props.interaction);
                if (props.debug) helper.debug.setProps(props.debug);
                if (props.handle) helper.handle.setProps(props.handle);

                if (cameraState.mode === 'orthographic') {
                    p.camera.stereo.name = 'off';
                }

                if (!doNotRequestDraw) {
                    requestDraw(true);
                }
            },
            getImagePass: (props: Partial<ImageProps> = {}) => {
                return new ImagePass(webgl, renderer, scene, camera, helper, passes.draw.wboitEnabled, props);
            },
            getRenderObjects(): GraphicsRenderObject[] {
                const renderObjects: GraphicsRenderObject[] = [];
                scene.forEach((_, ro) => renderObjects.push(ro));
                return renderObjects;
            },

            get props() {
                return getProps();
            },
            get input() {
                return input;
            },
            get stats() {
                return renderer.stats;
            },
            get interaction() {
                return interactionHelper.events;
            },
            dispose: () => {
                contextRestoredSub.unsubscribe();

                scene.clear();
                helper.debug.clear();
                controls.dispose();
                renderer.dispose();
                interactionHelper.dispose();
            }
        };

        function updateViewport() {
            const oldX = x, oldY = y, oldWidth = width, oldHeight = height;

            if (p.viewport.name === 'canvas') {
                x = 0;
                y = 0;
                width = gl.drawingBufferWidth;
                height = gl.drawingBufferHeight;
            } else if (p.viewport.name === 'static-frame') {
                x = p.viewport.params.x * webgl.pixelRatio;
                height = p.viewport.params.height * webgl.pixelRatio;
                y = gl.drawingBufferHeight - height - p.viewport.params.y * webgl.pixelRatio;
                width = p.viewport.params.width * webgl.pixelRatio;
            } else if (p.viewport.name === 'relative-frame') {
                x = Math.round(p.viewport.params.x * gl.drawingBufferWidth);
                height = Math.round(p.viewport.params.height * gl.drawingBufferHeight);
                y = Math.round(gl.drawingBufferHeight - height - p.viewport.params.y * gl.drawingBufferHeight);
                width = Math.round(p.viewport.params.width * gl.drawingBufferWidth);
            }

            if (oldX !== x || oldY !== y || oldWidth !== width || oldHeight !== height) {
                forceNextRender = true;
            }
        }

        function syncViewport() {
            pickHelper.setViewport(x, y, width, height);
            renderer.setViewport(x, y, width, height);
            Viewport.set(camera.viewport, x, y, width, height);
            Viewport.set(controls.viewport, x, y, width, height);
        }
    }
}