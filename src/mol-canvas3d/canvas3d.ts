/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Herman Bergwerf <post@hbergwerf.nl>
 */

import { BehaviorSubject, Subject, Subscription, debounceTime, merge } from 'rxjs';
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
import { AsyncPickData, DefaultPickOptions, PickData } from './passes/pick';
import { PickHelper } from './helper/pick-helper';
import { ImagePass, ImageProps } from './passes/image';
import { Sphere3D } from '../mol-math/geometry';
import { addConsoleStatsProvider, isDebugMode, isTimingMode, removeConsoleStatsProvider } from '../mol-util/debug';
import { CameraHelperParams } from './helper/camera-helper';
import { HandleHelperParams } from './helper/handle-helper';
import { StereoCamera, StereoCameraParams } from './camera/stereo';
import { Helper } from './helper/helper';
import { Passes } from './passes/passes';
import { shallowEqual } from '../mol-util';
import { MarkingParams } from './passes/marking';
import { degToRad, radToDeg } from '../mol-math/misc';
import { AssetManager } from '../mol-util/assets';
import { deepClone } from '../mol-util/object';
import { HiZParams, HiZPass } from './passes/hi-z';
import { IlluminationParams } from './passes/illumination';
import { isMobileBrowser } from '../mol-util/browser';
import { Ray3D } from '../mol-math/geometry/primitives/ray3d';
import { RayHelper } from './helper/ray-helper';
import { produce } from '../mol-util/produce';

export const CameraFogParams = {
    intensity: PD.Numeric(15, { min: 1, max: 100, step: 1 }),
};
export const Canvas3DParams = {
    camera: PD.Group({
        mode: PD.Select('perspective', PD.arrayToOptions(['perspective', 'orthographic'] as const), { label: 'Camera' }),
        helper: PD.Group(CameraHelperParams, { isFlat: true }),
        stereo: PD.MappedStatic('off', {
            on: PD.Group(StereoCameraParams),
            off: PD.Group({})
        }, { cycle: true, hideIf: p => p?.mode !== 'perspective' }),
        fov: PD.Numeric(45, { min: 10, max: 130, step: 1 }, { label: 'Field of View' }),
        scale: PD.Numeric(1, { min: 0.001, max: 1, step: 0.001 }, { label: 'Scene scale' }),
        manualReset: PD.Boolean(false, { isHidden: true }),
    }, { pivot: 'mode' }),
    cameraFog: PD.MappedStatic('on', {
        on: PD.Group(CameraFogParams),
        off: PD.Group({})
    }, { cycle: true, description: 'Show fog in the distance' }),
    cameraClipping: PD.Group({
        radius: PD.Numeric(100, { min: 0, max: 99, step: 1 }, { label: 'Clipping', description: 'How much of the scene to show.' }),
        far: PD.Boolean(true, { description: 'Hide scene in the distance' }),
        minNear: PD.Numeric(5, { min: 0.1, max: 100, step: 0.1 }, { description: 'Note, may cause performance issues rendering impostors when set too small and cause issues with outline rendering when too close to 0.' }),
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
    sceneRadiusFactor: PD.Numeric(1, { min: 1, max: 10, step: 0.1 }),
    transparentBackground: PD.Boolean(false),
    dpoitIterations: PD.Numeric(2, { min: 1, max: 10, step: 1 }),
    pickPadding: PD.Numeric(3, { min: 0, max: 10, step: 1 }, { description: 'Extra pixels to around target to check in case target is empty.' }),
    userInteractionReleaseMs: PD.Numeric(250, { min: 0, max: 1000, step: 1 }, { description: 'The time before the user is not considered interacting anymore.' }),

    multiSample: PD.Group(MultiSampleParams),
    postprocessing: PD.Group(PostprocessingParams),
    marking: PD.Group(MarkingParams),
    illumination: PD.Group(IlluminationParams),
    hiZ: PD.Group(HiZParams),
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
    readonly canvas?: HTMLCanvasElement
    readonly webgl: WebGLContext
    readonly input: InputObserver
    readonly passes: Passes
    readonly attribs: Readonly<Canvas3DContext.Attribs>
    readonly props: Readonly<Canvas3DContext.Props>
    readonly contextLost?: Subject<now.Timestamp>
    readonly contextRestored?: Subject<now.Timestamp>
    readonly assetManager: AssetManager
    readonly changed?: BehaviorSubject<undefined>
    readonly pixelScale: number

    syncPixelScale(): void
    setProps: (props?: Partial<Canvas3DContext.Props>) => void
    dispose: (options?: Partial<{ doNotForceWebGLContextLoss: boolean }>) => void
}

namespace Canvas3DContext {
    export const DefaultAttribs = {
        powerPreference: 'high-performance' as WebGLContextAttributes['powerPreference'],
        failIfMajorPerformanceCaveat: false,
        /** true by default to avoid issues with Safari (Jan 2021) */
        antialias: true,
        /** true to support multiple Canvas3D objects with a single context */
        preserveDrawingBuffer: true,
        preferWebGl1: false,

        handleResize: () => { },
    };
    export type Attribs = typeof DefaultAttribs

    export const Params = {
        resolutionMode: PD.Select('auto', PD.arrayToOptions(['auto', 'scaled', 'native'] as const)),
        pixelScale: PD.Numeric(1, { min: 0.1, max: 2, step: 0.05 }),
        pickScale: PD.Numeric(0.25, { min: 0.1, max: 1, step: 0.05 }),
        transparency: PD.Select('wboit', [['blended', 'Blended'], ['wboit', 'Weighted, Blended'], ['dpoit', 'Depth Peeling']] as const),
    };
    export const DefaultProps = PD.getDefaultValues(Params);
    export type Props = PD.Values<typeof Params>

    export function fromCanvas(canvas: HTMLCanvasElement, assetManager: AssetManager, attribs: Partial<Attribs> = {}, props: Partial<Props> = {}): Canvas3DContext {
        const a = { ...DefaultAttribs, ...attribs };
        const p = { ...DefaultProps, ...props };

        const { powerPreference, failIfMajorPerformanceCaveat, antialias, preserveDrawingBuffer, preferWebGl1 } = a;
        const gl = getGLContext(canvas, {
            powerPreference,
            failIfMajorPerformanceCaveat,
            antialias,
            preserveDrawingBuffer,
            alpha: true, // the renderer requires an alpha channel
            depth: true, // the renderer requires a depth buffer
            premultipliedAlpha: true, // the renderer outputs PMA
            preferWebGl1
        });
        if (gl === null) throw new Error('Could not create a WebGL rendering context');

        const getPixelScale = () => {
            const scaled = (p.pixelScale / (typeof window !== 'undefined' ? (window?.devicePixelRatio || 1) : 1));
            if (p.resolutionMode === 'auto') {
                return isMobileBrowser() ? scaled : p.pixelScale;
            }
            return p.resolutionMode === 'native' ? p.pixelScale : scaled;
        };
        const syncPixelScale = () => {
            const pixelScale = getPixelScale();
            input.setPixelScale(pixelScale);
            webgl.setPixelScale(pixelScale);
        };

        const { pickScale, transparency } = p;
        const pixelScale = getPixelScale();
        const input = InputObserver.fromElement(canvas, { pixelScale, preventGestures: true });
        const webgl = createContext(gl, { pixelScale });
        const passes = new Passes(webgl, assetManager, { pickScale, transparency });

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

        const contextLost = new Subject<now.Timestamp>();

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
                passes.pick.reset();
                passes.illumination.reset();
            });
            if (isDebugMode) console.log('context restored');
        };

        canvas.addEventListener('webglcontextlost', handleWebglContextLost, false);
        canvas.addEventListener('webglcontextrestored', handlewWebglContextRestored, false);

        const changed = new BehaviorSubject<undefined>(undefined);

        return {
            canvas,
            webgl,
            input,
            passes,
            attribs: a,
            get props() { return { ...p }; },
            contextLost,
            contextRestored: webgl.contextRestored,
            assetManager,
            changed,
            get pixelScale() { return getPixelScale(); },

            syncPixelScale,
            setProps: (props?: Partial<Props>) => {
                if (!props) return;

                let hasChanged = false;
                let pixelScaleNeedsUpdate = false;

                if (props.resolutionMode !== undefined && props.resolutionMode !== p.resolutionMode) {
                    p.resolutionMode = props.resolutionMode;
                    pixelScaleNeedsUpdate = true;
                }

                if (props.pixelScale !== undefined && props.pixelScale !== p.pixelScale) {
                    p.pixelScale = props.pixelScale;
                    pixelScaleNeedsUpdate = true;
                }

                if (pixelScaleNeedsUpdate) {
                    syncPixelScale();
                    a.handleResize();
                    hasChanged = true;
                }

                if (props.pickScale !== undefined && props.pickScale !== p.pickScale) {
                    p.pickScale = props.pickScale;
                    passes.setPickScale(props.pickScale);
                    hasChanged = true;
                }

                if (props.transparency !== undefined && props.transparency !== p.transparency) {
                    p.transparency = props.transparency;
                    passes.setTransparency(props.transparency);
                    hasChanged = true;
                }

                if (hasChanged) changed.next(undefined);
            },
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
    tick(t: now.Timestamp, options?: { isSynchronous?: boolean, manualDraw?: boolean, updateControls?: boolean }): void
    update(repr?: Representation.Any, keepBoundingSphere?: boolean): void
    clear(): void
    syncVisibility(): void

    requestDraw(): void

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
    identify(target: Vec2 | Ray3D): PickData | undefined
    asyncIdentify(target: Vec2 | Ray3D): AsyncPickData | undefined
    mark(loci: Representation.Loci, action: MarkerAction): void
    getLoci(pickingId: PickingId | undefined): Representation.Loci

    notifyDidDraw: boolean,
    readonly didDraw: BehaviorSubject<now.Timestamp>
    readonly commited: BehaviorSubject<now.Timestamp>
    readonly commitQueueSize: BehaviorSubject<number>
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

    export function create(ctx: Canvas3DContext, props: Partial<Canvas3DProps> = {}): Canvas3D {
        const { webgl, input, passes, assetManager, canvas, contextLost } = ctx;
        const p: Canvas3DProps = { ...deepClone(DefaultCanvas3DParams), ...deepClone(props) };

        const reprRenderObjects = new Map<Representation.Any, Set<GraphicsRenderObject>>();
        const reprUpdatedSubscriptions = new Map<Representation.Any, Subscription>();
        const reprCount = new BehaviorSubject(0);
        const interactionEvent = new Subject<void>();

        let startTime = now();
        const didDraw = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp);
        const commited = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp);
        const commitQueueSize = new BehaviorSubject<number>(0);

        const { contextRestored } = webgl;

        let x = 0;
        let y = 0;
        let width = 128;
        let height = 128;

        let forceNextRender = false;
        let currentTime = 0;

        updateViewport();
        const scene = Scene.create(webgl, passes.draw.transparency);

        function getSceneRadius() {
            return scene.boundingSphere.radius * p.sceneRadiusFactor;
        }

        const camera = new Camera({
            position: Vec3.create(0, 0, 100),
            mode: p.camera.mode,
            fog: p.cameraFog.name === 'on' ? p.cameraFog.params.intensity : 0,
            clipFar: p.cameraClipping.far,
            minNear: p.cameraClipping.minNear,
            fov: degToRad(p.camera.fov),
            scale: p.camera.scale,
        }, { x, y, width, height });
        const stereoCamera = new StereoCamera(camera, p.camera.stereo.params);

        const controls = TrackballControls.create(input, camera, scene, p.trackball);
        const helper = new Helper(webgl, scene, p);
        const hiZ = new HiZPass(webgl, passes.draw, canvas, p.hiZ);

        const renderer = Renderer.create(webgl, p.renderer);
        renderer.setOcclusionTest(hiZ.isOccluded);

        const pickOptions = {
            pickPadding: p.pickPadding,
            maxAsyncReadLag: DefaultPickOptions.maxAsyncReadLag,
        };
        const pickHelper = new PickHelper(webgl, renderer, scene, helper, passes.pick, { x, y, width, height }, pickOptions);
        const rayHelper = new RayHelper(webgl, renderer, scene, helper, pickOptions);
        const interactionHelper = new Canvas3dInteractionHelper(identify, asyncIdentify, getLoci, input, camera, controls, p.interaction);
        const multiSampleHelper = new MultiSampleHelper(passes.multiSample);

        passes.draw.postprocessing.background.update(camera, p.postprocessing.background, changed => {
            if (changed) requestDraw();
        });

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

        let markBuffer: [reprLoci: Representation.Loci, action: MarkerAction][] = [];

        function mark(reprLoci: Representation.Loci, action: MarkerAction) {
            // NOTE: might try to optimize a case with opposite actions for the
            //       same loci. Tho this might end up being more expensive (and error prone)
            //       then just applying everything "naively".
            markBuffer.push([reprLoci, action]);
        }

        function resolveMarking() {
            let changed = false;
            for (const [r, l] of markBuffer) {
                changed = applyMark(r, l) || changed;
            }
            markBuffer = [];
            if (changed) {
                scene.update(void 0, true);
                helper.handle.scene.update(void 0, true);
                helper.camera.scene.update(void 0, true);

                interactionEvent.next();
            }
            return changed;
        }

        function applyMark(reprLoci: Representation.Loci, action: MarkerAction) {
            const { repr, loci } = reprLoci;
            let changed = false;
            if (repr) {
                changed = repr.mark(loci, action) || changed;
            } else {
                reprRenderObjects.forEach((_, _repr) => { changed = _repr.mark(loci, action) || changed; });
            }
            changed = helper.handle.mark(loci, action) || changed;
            changed = helper.camera.mark(loci, action) || changed;
            return changed;
        }

        let fenceSync: WebGLSync | null = null;

        function render(force: boolean) {
            if (webgl.isContextLost) return false;

            let resized = false;
            if (resizeRequested) {
                handleResize(false);
                resizeRequested = false;
                resized = true;
            }

            const drs = webgl.getDrawingBufferSize();
            if (x > drs.width || x + width < 0 ||
                y > drs.height || y + height < 0
            ) return false;

            if (fenceSync !== null) {
                if (webgl.checkSyncStatus(fenceSync)) {
                    fenceSync = null;
                } else {
                    return false;
                }
            }

            const markingUpdated = resolveMarking() && (renderer.props.colorMarker || p.marking.enabled);

            let didRender = false;
            controls.update(currentTime);
            const cameraChanged = camera.update();

            const shouldRender = force || cameraChanged || resized || forceNextRender;
            forceNextRender = false;

            if (passes.illumination.supported && p.illumination.enabled) {
                if (shouldRender || markingUpdated) {
                    renderer.setOcclusionTest(null);
                    passes.illumination.restart();
                }

                if (passes.illumination.shouldRender(p)
                    && ((!isActivelyInteracting && scene.count > 0) || passes.illumination.iteration === 0 || p.userInteractionReleaseMs === 0)
                ) {
                    if (isTimingMode) webgl.timer.mark('Canvas3D.render', { captureStats: true });
                    const ctx = { renderer, camera, scene, helper };
                    passes.illumination.render(ctx, p, true);
                    if (isTimingMode) webgl.timer.markEnd('Canvas3D.render');

                    // if only marking has updated, do not set the flag to dirty
                    pickHelper.dirty = pickHelper.dirty || shouldRender;
                    didRender = true;
                }
            } else {
                const multiSampleChanged = multiSampleHelper.update(markingUpdated || shouldRender, p.multiSample);

                if (shouldRender || multiSampleChanged || markingUpdated) {
                    renderer.setOcclusionTest(hiZ.isOccluded);

                    let cam: Camera | StereoCamera = camera;
                    if (p.camera.stereo.name === 'on') {
                        stereoCamera.update();
                        cam = stereoCamera;
                    }

                    if (isTimingMode) webgl.timer.mark('Canvas3D.render', { captureStats: true });
                    const ctx = { renderer, camera: cam, scene, helper };
                    if (MultiSamplePass.isEnabled(p.multiSample)) {
                        const forceOn = p.multiSample.reduceFlicker && !cameraChanged && markingUpdated && !controls.isAnimating;
                        multiSampleHelper.render(ctx, p, true, forceOn);
                    } else {
                        passes.draw.render(ctx, p, true);
                    }
                    hiZ.render(camera);
                    if (isTimingMode) webgl.timer.markEnd('Canvas3D.render');

                    // if only marking has updated, do not set the flag to dirty
                    pickHelper.dirty = pickHelper.dirty || shouldRender;
                    didRender = true;
                }
            }

            if (didRender) {
                fenceSync = webgl.getFenceSync();
            }

            return didRender;
        }

        let forceDrawAfterAllCommited = false;
        let drawPaused = false;
        let isContextLost = false;

        function draw(options?: { force?: boolean }) {
            if (drawPaused || isContextLost) return;
            if (render(!!options?.force) && notifyDidDraw) {
                didDraw.next(now() - startTime as now.Timestamp);
            }
        }

        function requestDraw() {
            forceNextRender = true;
        }

        let animationFrameHandle = 0;

        function tick(t: now.Timestamp, options?: { isSynchronous?: boolean, manualDraw?: boolean, updateControls?: boolean }) {
            if (isContextLost) return;

            currentTime = t;
            commit(options?.isSynchronous);

            // update the controler before the camera transition
            if (options?.updateControls) {
                controls.update(currentTime);
            }

            camera.transition.tick(currentTime);
            hiZ.tick();

            if (options?.manualDraw) {
                return;
            }

            draw();
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

        function identify(target: Vec2 | Ray3D): PickData | undefined {
            if (webgl.isContextLost) return undefined;

            if ('origin' in target) {
                return rayHelper.identify(target, camera);
            } else {
                const cam = (p.camera.stereo.name === 'on') ? stereoCamera : camera;
                return pickHelper.identify(target[0], target[1], cam);
            }
        }

        function asyncIdentify(target: Vec2 | Ray3D): AsyncPickData | undefined {
            if (webgl.isContextLost) return undefined;

            if ('origin' in target) {
                return rayHelper.asyncIdentify(target, camera);
            } else {
                const cam = (p.camera.stereo.name === 'on') ? stereoCamera : camera;
                return pickHelper.asyncIdentify(target[0], target[1], cam);
            }
        }

        function commit(isSynchronous: boolean = false) {
            const allCommited = commitScene(isSynchronous);
            // Only reset the camera after the full scene has been commited.
            if (allCommited) {
                resolveCameraReset();
                if (forceDrawAfterAllCommited) {
                    if (helper.debug.isEnabled) helper.debug.update();
                    draw({ force: true });
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
                camera.setState({ ...snapshot, radiusMax: getSceneRadius() }, duration);
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

            let cameraSphereOverlapsNone = true, isEmpty = true;
            Sphere3D.set(cameraSphere, camera.state.target, camera.state.radius);

            // check if any renderable has moved outside of the old bounding sphere
            // and if no renderable is overlapping with the camera sphere
            for (const r of scene.renderables) {
                if (!r.state.visible) continue;

                const b = r.values.boundingSphere.ref.value;
                if (!b.radius) continue;

                isEmpty = false;
                const cameraDist = Vec3.distance(cameraSphere.center, b.center);
                if ((cameraDist > cameraSphere.radius || cameraDist > b.radius || b.radius > camera.state.radiusMax) && !Sphere3D.includes(oldBoundingSphereVisible, b)) return true;
                if (Sphere3D.overlaps(cameraSphere, b)) cameraSphereOverlapsNone = false;
            }

            return cameraSphereOverlapsNone || (!isEmpty && cameraSphere.radius <= 0.1);
        }

        const sceneCommitTimeoutMs = 250;
        function commitScene(isSynchronous: boolean) {
            if (!scene.needsCommit) return true;

            // snapshot the current bounding sphere of visible objects
            Sphere3D.copy(oldBoundingSphereVisible, scene.boundingSphereVisible);

            // clear hi-Z buffer when scene changes
            hiZ.clear();

            if (!scene.commit(isSynchronous ? void 0 : sceneCommitTimeoutMs)) {
                commitQueueSize.next(scene.commitQueueSize);
                return false;
            }
            commitQueueSize.next(0);

            if (helper.debug.isEnabled) helper.debug.update();
            if (!p.camera.manualReset && (reprCount.value === 0 || shouldResetCamera())) {
                cameraResetRequested = true;
            }
            if (oldBoundingSphereVisible.radius === 0) nextCameraResetDuration = 0;

            if (!p.camera.manualReset) camera.setState({ radiusMax: getSceneRadius() }, 0);
            reprCount.next(reprRenderObjects.size);
            if (isDebugMode) consoleStats();

            return true;
        }

        function consoleStats() {
            const items = scene.renderables.map(r => ({
                drawCount: r.values.drawCount.ref.value,
                instanceCount: r.values.instanceCount.ref.value,
                materialId: r.materialId,
                renderItemId: r.id,
            }));

            console.groupCollapsed(`${items.length} RenderItems`);

            if (items.length < 50) {
                console.table(items);
            } else {
                console.log(items);
            }
            console.log(JSON.stringify(webgl.stats, undefined, 4));

            const { texture, attribute, elements } = webgl.resources.getByteCounts();
            console.log(JSON.stringify({
                texture: `${(texture / 1024 / 1024).toFixed(3)} MiB`,
                attribute: `${(attribute / 1024 / 1024).toFixed(3)} MiB`,
                elements: `${(elements / 1024 / 1024).toFixed(3)} MiB`,
            }, undefined, 4));

            console.log(JSON.stringify(webgl.timer.formatedStats(), undefined, 4));

            console.groupEnd();
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
                ? 100 - Math.round((camera.transition.target.radius / getSceneRadius()) * 100)
                : 0;

            return {
                camera: {
                    mode: camera.state.mode,
                    helper: { ...helper.camera.props },
                    stereo: { ...p.camera.stereo },
                    fov: Math.round(radToDeg(camera.state.fov)),
                    scale: camera.state.scale,
                    manualReset: !!p.camera.manualReset
                },
                cameraFog: camera.state.fog > 0
                    ? { name: 'on' as const, params: { intensity: camera.state.fog } }
                    : { name: 'off' as const, params: {} },
                cameraClipping: { far: camera.state.clipFar, radius, minNear: camera.state.minNear },
                cameraResetDurationMs: p.cameraResetDurationMs,
                sceneRadiusFactor: p.sceneRadiusFactor,
                transparentBackground: p.transparentBackground,
                dpoitIterations: p.dpoitIterations,
                pickPadding: p.pickPadding,
                userInteractionReleaseMs: p.userInteractionReleaseMs,
                viewport: p.viewport,

                postprocessing: { ...p.postprocessing },
                marking: { ...p.marking },
                multiSample: { ...p.multiSample },
                illumination: { ...p.illumination },
                hiZ: { ...hiZ.props },
                renderer: { ...renderer.props },
                trackball: { ...controls.props },
                interaction: { ...interactionHelper.props },
                debug: { ...helper.debug.props },
                handle: { ...helper.handle.props },
            };
        }

        const contextLostSub = contextLost?.subscribe(() => {
            isContextLost = true;
            fenceSync = null;
            pickHelper.dirty = true;
        });

        const contextRestoredSub = contextRestored.subscribe(() => {
            pickHelper.reset();
            rayHelper.reset();
            hiZ.reset();

            scene.forEach(r => {
                if (r.values.meta?.ref.value.reset) {
                    r.values.meta.ref.value.reset();
                    r.update();
                }
            });

            isContextLost = false;

            draw({ force: true });
            // Unclear why, but in Chrome with wboit enabled the first `draw` only clears
            // the drawingBuffer. Note that in Firefox the drawingBuffer is preserved after
            // context loss so it is unclear if it behaves the same.
            draw({ force: true });
        });

        const resized = new BehaviorSubject<any>(0);

        function handleResize(draw = true) {
            passes.updateSize();
            updateViewport();
            syncViewport();
            if (draw) requestDraw();
            resized.next(+new Date());
        }

        addConsoleStatsProvider(consoleStats);

        const ctxChangedSub = ctx.changed?.subscribe(() => {
            scene.setTransparency(passes.draw.transparency);
            requestDraw();
        });

        // Monitor user interactions
        let isDragging = false;
        let isActivelyInteracting = false;
        let interactionSubs = [
            input.drag.subscribe(() => {
                isDragging = true;
            }),
            input.interactionEnd.subscribe(() => {
                isDragging = false;
            }),
            merge(
                input.drag,
                input.pinch,
                input.wheel,
                input.interactionEnd,
            ).subscribe(() => {
                interactionEvent.next();
            }),
            interactionEvent.subscribe(() => {
                isActivelyInteracting = true;
            }),
            interactionEvent.pipe(
                debounceTime(p.userInteractionReleaseMs)
            ).subscribe(() => {
                isActivelyInteracting = isDragging;
                if (!isDragging && passes.illumination.supported && p.illumination.enabled) {
                    requestDraw();
                }
            }),
        ];

        //

        if (isDebugMode && canvas) {
            let occlusionLoci: Loci | undefined = undefined;

            const printOcclusion = (loci: Loci | undefined) => {
                const s = loci && Loci.getBoundingSphere(Loci.normalize(loci, 'residue'));
                hiZ.debugOcclusion(s);
            };

            input.click.subscribe(e => {
                if (!e.modifiers.control || e.button !== 2) return;

                const p = identify(Vec2.create(e.x, e.y));
                if (!p) {
                    occlusionLoci = undefined;
                    printOcclusion(occlusionLoci);
                    return;
                }

                const l = getLoci(p.id);
                occlusionLoci = l.loci;
                printOcclusion(occlusionLoci);
            });

            didDraw.subscribe(() => {
                setTimeout(() => {
                    printOcclusion(occlusionLoci);
                }, 100);
            });
        }

        //

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
                if (fenceSync !== null) {
                    webgl.deleteSync(fenceSync);
                    fenceSync = null;
                }
                requestDraw();
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
                requestDraw();
            },

            requestDraw,
            tick,
            animate,
            resetTime,
            pause,
            resume: () => { drawPaused = false; },
            identify,
            asyncIdentify,
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
            commitQueueSize,
            reprCount,
            resized,
            setProps: (properties, doNotRequestDraw = false) => {
                const props: PartialCanvas3DProps = typeof properties === 'function'
                    ? produce(getProps(), properties as any)
                    : properties;

                if (props.sceneRadiusFactor !== undefined) {
                    p.sceneRadiusFactor = props.sceneRadiusFactor;
                    camera.setState({ radiusMax: getSceneRadius() }, 0);
                }

                const cameraState: Partial<Camera.Snapshot> = Object.create(null);
                if (props.camera && props.camera.mode !== undefined && props.camera.mode !== camera.state.mode) {
                    cameraState.mode = props.camera.mode;
                }
                const oldFov = Math.round(radToDeg(camera.state.fov));
                if (props.camera && props.camera.fov !== undefined && props.camera.fov !== oldFov) {
                    cameraState.fov = degToRad(props.camera.fov);
                }
                if (props.camera && props.camera.scale !== undefined && props.camera.scale !== cameraState.scale) {
                    cameraState.scale = props.camera.scale;
                }
                if (props.cameraFog !== undefined && props.cameraFog.params) {
                    const newFog = props.cameraFog.name === 'on' ? props.cameraFog.params.intensity : 0;
                    if (newFog !== camera.state.fog) cameraState.fog = newFog;
                }
                if (props.cameraClipping !== undefined) {
                    if (props.cameraClipping.far !== undefined && props.cameraClipping.far !== camera.state.clipFar) {
                        cameraState.clipFar = props.cameraClipping.far;
                    }
                    if (props.cameraClipping.minNear !== undefined && props.cameraClipping.minNear !== camera.state.minNear) {
                        cameraState.minNear = props.cameraClipping.minNear;
                    }
                    if (props.cameraClipping.radius !== undefined) {
                        const radius = (getSceneRadius() / 100) * (100 - props.cameraClipping.radius);
                        if (radius > 0 && radius !== cameraState.radius) {
                            // if radius = 0, NaNs happen
                            cameraState.radius = Math.max(radius, 0.01);
                        }
                    }
                }
                if (Object.keys(cameraState).length > 0) camera.setState(cameraState);

                if (props.camera?.helper) helper.camera.setProps(props.camera.helper);
                if (props.camera?.manualReset !== undefined) p.camera.manualReset = props.camera.manualReset;
                if (props.camera?.stereo !== undefined) {
                    Object.assign(p.camera.stereo, props.camera.stereo);
                    stereoCamera.setProps(p.camera.stereo.params);
                }
                if (props.cameraResetDurationMs !== undefined) p.cameraResetDurationMs = props.cameraResetDurationMs;
                if (props.transparentBackground !== undefined) p.transparentBackground = props.transparentBackground;
                if (props.dpoitIterations !== undefined) p.dpoitIterations = props.dpoitIterations;
                if (props.pickPadding !== undefined) {
                    p.pickPadding = props.pickPadding;
                    pickHelper.setPickPadding(p.pickPadding);
                }
                if (props.userInteractionReleaseMs !== undefined) p.userInteractionReleaseMs = props.userInteractionReleaseMs;
                if (props.viewport !== undefined) {
                    const doNotUpdate = p.viewport === props.viewport ||
                        (p.viewport.name === props.viewport.name && shallowEqual(p.viewport.params, props.viewport.params));

                    if (!doNotUpdate) {
                        p.viewport = props.viewport;
                        updateViewport();
                        syncViewport();
                    }
                }

                if (props.postprocessing?.background) {
                    Object.assign(p.postprocessing.background, props.postprocessing.background);
                    passes.draw.postprocessing.background.update(camera, p.postprocessing.background, changed => {
                        if (changed && !doNotRequestDraw) requestDraw();
                    });
                }
                if (props.postprocessing) Object.assign(p.postprocessing, props.postprocessing);
                if (props.marking) Object.assign(p.marking, props.marking);
                if (props.illumination) Object.assign(p.illumination, props.illumination);
                if (props.multiSample) Object.assign(p.multiSample, props.multiSample);
                if (props.hiZ) hiZ.setProps(props.hiZ);
                if (props.renderer) renderer.setProps(props.renderer);
                if (props.trackball) controls.setProps(props.trackball);
                if (props.interaction) interactionHelper.setProps(props.interaction);
                if (props.debug) helper.debug.setProps(props.debug);
                if (props.handle) helper.handle.setProps(props.handle);

                if (cameraState.mode === 'orthographic') {
                    p.camera.stereo.name = 'off';
                }

                if (!doNotRequestDraw) {
                    requestDraw();
                }
            },
            getImagePass: (props: Partial<ImageProps> = {}) => {
                return new ImagePass(webgl, assetManager, renderer, scene, camera, helper, passes.draw.transparency, props);
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
                contextLostSub?.unsubscribe();
                contextRestoredSub.unsubscribe();
                ctxChangedSub?.unsubscribe();

                for (const s of interactionSubs) s.unsubscribe();
                interactionSubs = [];

                cancelAnimationFrame(animationFrameHandle);

                markBuffer = [];

                scene.clear();
                helper.debug.clear();
                controls.dispose();
                renderer.dispose();
                interactionHelper.dispose();
                hiZ.dispose();
                pickHelper.dispose();
                rayHelper.dispose();

                if (fenceSync !== null) {
                    webgl.deleteSync(fenceSync);
                    fenceSync = null;
                }

                removeConsoleStatsProvider(consoleStats);
            }
        };

        function updateViewport() {
            const oldX = x, oldY = y, oldWidth = width, oldHeight = height;
            const drs = webgl.getDrawingBufferSize();

            if (p.viewport.name === 'canvas') {
                x = 0;
                y = 0;
                width = drs.width;
                height = drs.height;
            } else if (p.viewport.name === 'static-frame') {
                x = p.viewport.params.x * webgl.pixelRatio;
                height = p.viewport.params.height * webgl.pixelRatio;
                y = drs.height - height - p.viewport.params.y * webgl.pixelRatio;
                width = p.viewport.params.width * webgl.pixelRatio;
            } else if (p.viewport.name === 'relative-frame') {
                x = Math.round(p.viewport.params.x * drs.width);
                height = Math.round(p.viewport.params.height * drs.height);
                y = Math.round(drs.height - height - p.viewport.params.y * drs.height);
                width = Math.round(p.viewport.params.width * drs.width);
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
            hiZ.setViewport(x, y, width, height);
        }
    }
}
