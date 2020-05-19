/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, Subscription } from 'rxjs';
import { now } from '../mol-util/now';
import { Vec3 } from '../mol-math/linear-algebra';
import InputObserver, { ModifiersKeys, ButtonsType } from '../mol-util/input/input-observer';
import Renderer, { RendererStats, RendererParams } from '../mol-gl/renderer';
import { GraphicsRenderObject } from '../mol-gl/render-object';
import { TrackballControls, TrackballControlsParams } from './controls/trackball';
import { Viewport } from './camera/util';
import { createContext, WebGLContext, getGLContext } from '../mol-gl/webgl/context';
import { Representation } from '../mol-repr/representation';
import Scene from '../mol-gl/scene';
import { GraphicsRenderVariant } from '../mol-gl/webgl/render-item';
import { PickingId } from '../mol-geo/geometry/picking';
import { MarkerAction } from '../mol-util/marker-action';
import { Loci, EmptyLoci, isEmptyLoci } from '../mol-model/loci';
import { Camera } from './camera';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { BoundingSphereHelper, DebugHelperParams } from './helper/bounding-sphere-helper';
import { SetUtils } from '../mol-util/set';
import { Canvas3dInteractionHelper } from './helper/interaction-events';
import { PostprocessingParams, PostprocessingPass } from './passes/postprocessing';
import { MultiSampleParams, MultiSamplePass } from './passes/multi-sample';
import { PixelData } from '../mol-util/image';
import { readTexture } from '../mol-gl/compute/util';
import { DrawPass } from './passes/draw';
import { PickPass } from './passes/pick';
import { ImagePass, ImageProps } from './passes/image';
import { Sphere3D } from '../mol-math/geometry';
import { isDebugMode } from '../mol-util/debug';
import { CameraHelperParams } from './helper/camera-helper';
import { produce } from 'immer';
import { HandleHelper, HandleHelperParams } from './helper/handle-helper';

export const Canvas3DParams = {
    camera: PD.Group({
        mode: PD.Select('perspective', [['perspective', 'Perspective'], ['orthographic', 'Orthographic']] as const, { label: 'Camera' }),
        helper: PD.Group(CameraHelperParams, { isFlat: true })
    }, { pivot: 'mode' }),
    cameraFog: PD.MappedStatic('on', {
        on: PD.Group({
            intensity: PD.Numeric(50, { min: 1, max: 100, step: 1 }),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Show fog in the distance' }),
    cameraClipping: PD.Group({
        radius: PD.Numeric(100, { min: 0, max: 99, step: 1 }, { label: 'Clipping', description: 'How much of the scene to show.' }),
        far: PD.Boolean(true, { description: 'Hide scene in the distance' }),
    }, { pivot: 'radius' }),

    cameraResetDurationMs: PD.Numeric(250, { min: 0, max: 1000, step: 1 }, { description: 'The time it takes to reset the camera.' }),
    transparentBackground: PD.Boolean(false),

    multiSample: PD.Group(MultiSampleParams),
    postprocessing: PD.Group(PostprocessingParams),
    renderer: PD.Group(RendererParams),
    trackball: PD.Group(TrackballControlsParams),
    debug: PD.Group(DebugHelperParams),
    handle: PD.Group(HandleHelperParams),
};
export const DefaultCanvas3DParams = PD.getDefaultValues(Canvas3DParams);
export type Canvas3DProps = PD.Values<typeof Canvas3DParams>

export { Canvas3D };

interface Canvas3D {
    readonly webgl: WebGLContext,

    add(repr: Representation.Any): void
    remove(repr: Representation.Any): void
    /**
     * This function must be called if animate() is not set up so that add/remove actions take place.
     */
    commit(isSynchronous?: boolean): void
    update(repr?: Representation.Any, keepBoundingSphere?: boolean): void
    clear(): void
    syncVisibility(): void

    requestDraw(force?: boolean): void
    animate(): void
    identify(x: number, y: number): PickingId | undefined
    mark(loci: Representation.Loci, action: MarkerAction): void
    getLoci(pickingId: PickingId): Representation.Loci

    readonly didDraw: BehaviorSubject<now.Timestamp>
    readonly reprCount: BehaviorSubject<number>

    handleResize(): void
    /** Focuses camera on scene's bounding sphere, centered and zoomed. */
    requestCameraReset(options?: { durationMs?: number, snapshot?: Partial<Camera.Snapshot> }): void
    readonly camera: Camera
    readonly boundingSphere: Readonly<Sphere3D>
    getPixelData(variant: GraphicsRenderVariant): PixelData
    setProps(props: Partial<Canvas3DProps> | ((old: Canvas3DProps) => Partial<Canvas3DProps> | void)): void
    getImagePass(props: Partial<ImageProps>): ImagePass

    /** Returns a copy of the current Canvas3D instance props */
    readonly props: Readonly<Canvas3DProps>
    readonly input: InputObserver
    readonly stats: RendererStats
    readonly interaction: Canvas3dInteractionHelper['events']

    dispose(): void
}

const requestAnimationFrame = typeof window !== 'undefined' ? window.requestAnimationFrame : (f: (time: number) => void) => setImmediate(()=>f(Date.now()));

namespace Canvas3D {
    export interface HoverEvent { current: Representation.Loci, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys }
    export interface ClickEvent { current: Representation.Loci, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys }

    export function fromCanvas(canvas: HTMLCanvasElement, props: Partial<Canvas3DProps> = {}) {
        const gl = getGLContext(canvas, {
            alpha: true,
            antialias: true,
            depth: true,
            preserveDrawingBuffer: true,
            premultipliedAlpha: false,
        });
        if (gl === null) throw new Error('Could not create a WebGL rendering context');
        const input = InputObserver.fromElement(canvas);
        const webgl = createContext(gl);

        if (isDebugMode) {
            const loseContextExt = gl.getExtension('WEBGL_lose_context');
            if (loseContextExt) {
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

        canvas.addEventListener('webglcontextlost', e => {
            webgl.setContextLost();
            e.preventDefault();
            if (isDebugMode) console.log('context lost');
        }, false);

        canvas.addEventListener('webglcontextrestored', () => {
            if (!webgl.isContextLost) return;
            webgl.handleContextRestored();
            if (isDebugMode) console.log('context restored');
        }, false);

        return Canvas3D.create(webgl, input, props);
    }

    export function create(webgl: WebGLContext, input: InputObserver, props: Partial<Canvas3DProps> = {}): Canvas3D {
        const p = { ...DefaultCanvas3DParams, ...props };

        const reprRenderObjects = new Map<Representation.Any, Set<GraphicsRenderObject>>();
        const reprUpdatedSubscriptions = new Map<Representation.Any, Subscription>();
        const reprCount = new BehaviorSubject(0);

        const startTime = now();
        const didDraw = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp);

        const { gl, contextRestored } = webgl;

        let width = gl.drawingBufferWidth;
        let height = gl.drawingBufferHeight;

        const scene = Scene.create(webgl);

        const camera = new Camera({
            position: Vec3.create(0, 0, 100),
            mode: p.camera.mode,
            fog: p.cameraFog.name === 'on' ? p.cameraFog.params.intensity : 0,
            clipFar: p.cameraClipping.far
        });

        const controls = TrackballControls.create(input, camera, p.trackball);
        const renderer = Renderer.create(webgl, p.renderer);
        const debugHelper = new BoundingSphereHelper(webgl, scene, p.debug);
        const handleHelper = new HandleHelper(webgl, p.handle);
        const interactionHelper = new Canvas3dInteractionHelper(identify, getLoci, input);

        const drawPass = new DrawPass(webgl, renderer, scene, camera, debugHelper, handleHelper, {
            cameraHelper: p.camera.helper
        });
        const pickPass = new PickPass(webgl, renderer, scene, camera, handleHelper, 0.5);
        const postprocessing = new PostprocessingPass(webgl, camera, drawPass, p.postprocessing);
        const multiSample = new MultiSamplePass(webgl, camera, drawPass, postprocessing, p.multiSample);

        let drawPending = false;
        let cameraResetRequested = false;
        let nextCameraResetDuration: number | undefined = void 0;
        let nextCameraResetSnapshot: Partial<Camera.Snapshot> | undefined = void 0;

        function getLoci(pickingId: PickingId) {
            let loci: Loci = EmptyLoci;
            let repr: Representation.Any = Representation.Empty;
            loci = handleHelper.getLoci(pickingId);
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
            return { loci, repr };
        }

        function mark(reprLoci: Representation.Loci, action: MarkerAction) {
            const { repr, loci } = reprLoci;
            let changed = false;
            if (repr) {
                changed = repr.mark(loci, action);
            } else {
                changed = handleHelper.mark(loci, action);
                reprRenderObjects.forEach((_, _repr) => { changed = _repr.mark(loci, action) || changed; });
            }
            if (changed) {
                scene.update(void 0, true);
                handleHelper.scene.update(void 0, true);
                const prevPickDirty = pickPass.pickDirty;
                draw(true);
                pickPass.pickDirty = prevPickDirty; // marking does not change picking buffers
            }
        }

        function render(force: boolean) {
            if (webgl.isContextLost) return false;

            let didRender = false;
            controls.update(currentTime);
            Viewport.set(camera.viewport, 0, 0, width, height);
            const cameraChanged = camera.update();
            multiSample.update(force || cameraChanged, currentTime);

            if (force || cameraChanged || multiSample.enabled) {
                renderer.setViewport(0, 0, width, height);
                if (multiSample.enabled) {
                    multiSample.render(true, p.transparentBackground);
                } else {
                    drawPass.render(!postprocessing.enabled, p.transparentBackground);
                    if (postprocessing.enabled) postprocessing.render(true);
                }
                pickPass.pickDirty = true;
                didRender = true;
            }

            return didRender;
        }

        let forceNextDraw = false;
        let forceDrawAfterAllCommited = false;
        let currentTime = 0;

        function draw(force?: boolean) {
            if (render(!!force || forceNextDraw)) {
                didDraw.next(now() - startTime as now.Timestamp);
            }
            forceNextDraw = false;
            drawPending = false;
        }

        function requestDraw(force?: boolean) {
            if (drawPending) return;
            drawPending = true;
            forceNextDraw = !!force;
        }

        function animate() {
            currentTime = now();
            commit();
            camera.transition.tick(currentTime);

            draw(false);
            if (!camera.transition.inTransition && !webgl.isContextLost) {
                interactionHelper.tick(currentTime);
            }
            requestAnimationFrame(animate);
        }

        function identify(x: number, y: number): PickingId | undefined {
            return webgl.isContextLost ? undefined : pickPass.identify(x, y);
        }

        function commit(isSynchronous: boolean = false) {
            const allCommited = commitScene(isSynchronous);
            // Only reset the camera after the full scene has been commited.
            if (allCommited) {
                resolveCameraReset();
                if (forceDrawAfterAllCommited) {
                    draw(true);
                    forceDrawAfterAllCommited = false;
                }
            }
        }

        function resolveCameraReset() {
            if (!cameraResetRequested) return;

            const { center, radius } = scene.boundingSphereVisible;
            if (radius > 0) {
                const duration = nextCameraResetDuration === undefined ? p.cameraResetDurationMs : nextCameraResetDuration;
                const focus = camera.getFocus(center, radius);
                const snapshot = nextCameraResetSnapshot ? { ...focus, ...nextCameraResetSnapshot } : focus;
                camera.setState(snapshot, duration);
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

                if (!Sphere3D.includes(oldBoundingSphereVisible, b)) return true;
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

            if (debugHelper.isEnabled) debugHelper.update();
            if (reprCount.value === 0 || shouldResetCamera()) {
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
            })));
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
                scene.update(repr.renderObjects, false, true);
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
                    helper: { ...drawPass.props.cameraHelper }
                },
                cameraFog: camera.state.fog > 0
                    ? { name: 'on' as const, params: { intensity: camera.state.fog } }
                    : { name: 'off' as const, params: {} },
                cameraClipping: { far: camera.state.clipFar, radius },
                cameraResetDurationMs: p.cameraResetDurationMs,
                transparentBackground: p.transparentBackground,

                postprocessing: { ...postprocessing.props },
                multiSample: { ...multiSample.props },
                renderer: { ...renderer.props },
                trackball: { ...controls.props },
                debug: { ...debugHelper.props },
                handle: { ...handleHelper.props },
            };
        }

        handleResize();

        const contextRestoredSub = contextRestored.subscribe(() => {
            pickPass.pickDirty = true;
            draw(true);
        });

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
                debugHelper.clear();
                requestDraw(true);
                reprCount.next(reprRenderObjects.size);
            },
            syncVisibility: () => {
                if (camera.state.radiusMax === 0) {
                    cameraResetRequested = true;
                    nextCameraResetDuration = 0;
                }

                if (scene.syncVisibility()) {
                    if (debugHelper.isEnabled) debugHelper.update();
                }
                requestDraw(true);
            },

            // draw,
            requestDraw,
            animate,
            identify,
            mark,
            getLoci,

            handleResize,
            requestCameraReset: options => {
                nextCameraResetDuration = options?.durationMs;
                nextCameraResetSnapshot = options?.snapshot;
                cameraResetRequested = true;
            },
            camera,
            boundingSphere: scene.boundingSphere,
            getPixelData: (variant: GraphicsRenderVariant) => {
                switch (variant) {
                    case 'color': return webgl.getDrawingBufferPixelData();
                    case 'pickObject': return pickPass.objectPickTarget.getPixelData();
                    case 'pickInstance': return pickPass.instancePickTarget.getPixelData();
                    case 'pickGroup': return pickPass.groupPickTarget.getPixelData();
                    case 'depth': return readTexture(webgl, drawPass.depthTexture) as PixelData;
                }
            },
            didDraw,
            reprCount,
            setProps: (properties) => {
                const props: Partial<Canvas3DProps> = typeof properties === 'function'
                    ? produce(getProps(), properties)
                    : properties;

                const cameraState: Partial<Camera.Snapshot> = Object.create(null);
                if (props.camera && props.camera.mode !== undefined && props.camera.mode !== camera.state.mode) {
                    cameraState.mode = props.camera.mode;
                }
                if (props.cameraFog !== undefined) {
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

                if (props.camera?.helper) drawPass.setProps({ cameraHelper: props.camera.helper });
                if (props.cameraResetDurationMs !== undefined) p.cameraResetDurationMs = props.cameraResetDurationMs;
                if (props.transparentBackground !== undefined) p.transparentBackground = props.transparentBackground;

                if (props.postprocessing) postprocessing.setProps(props.postprocessing);
                if (props.multiSample) multiSample.setProps(props.multiSample);
                if (props.renderer) renderer.setProps(props.renderer);
                if (props.trackball) controls.setProps(props.trackball);
                if (props.debug) debugHelper.setProps(props.debug);
                if (props.handle) handleHelper.setProps(props.handle);

                requestDraw(true);
            },
            getImagePass: (props: Partial<ImageProps> = {}) => {
                return new ImagePass(webgl, renderer, scene, camera, debugHelper, handleHelper, props);
            },

            get props() {
                const radius = scene.boundingSphere.radius > 0
                    ? 100 - Math.round((camera.transition.target.radius / scene.boundingSphere.radius) * 100)
                    : 0;

                return {
                    camera: {
                        mode: camera.state.mode,
                        helper: { ...drawPass.props.cameraHelper }
                    },
                    cameraFog: camera.state.fog > 0
                        ? { name: 'on' as const, params: { intensity: camera.state.fog } }
                        : { name: 'off' as const, params: {} },
                    cameraClipping: { far: camera.state.clipFar, radius },
                    cameraResetDurationMs: p.cameraResetDurationMs,
                    transparentBackground: p.transparentBackground,

                    postprocessing: { ...postprocessing.props },
                    multiSample: { ...multiSample.props },
                    renderer: { ...renderer.props },
                    trackball: { ...controls.props },
                    debug: { ...debugHelper.props },
                    handle: { ...handleHelper.props },
                };
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
                debugHelper.clear();
                input.dispose();
                controls.dispose();
                renderer.dispose();
                interactionHelper.dispose();
            }
        };

        function handleResize() {
            width = gl.drawingBufferWidth;
            height = gl.drawingBufferHeight;

            renderer.setViewport(0, 0, width, height);
            Viewport.set(camera.viewport, 0, 0, width, height);
            Viewport.set(controls.viewport, 0, 0, width, height);

            drawPass.setSize(width, height);
            pickPass.setSize(width, height);
            postprocessing.setSize(width, height);
            multiSample.setSize(width, height);

            requestDraw(true);
        }
    }
}