/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BehaviorSubject, Subscription } from 'rxjs';
import { now } from '../mol-util/now';
import { Vec3 } from '../mol-math/linear-algebra'
import InputObserver, { ModifiersKeys, ButtonsType } from '../mol-util/input/input-observer'
import Renderer, { RendererStats, RendererParams } from '../mol-gl/renderer'
import { GraphicsRenderObject } from '../mol-gl/render-object'
import { TrackballControls, TrackballControlsParams } from './controls/trackball'
import { Viewport } from './camera/util'
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
import { GLRenderingContext } from '../mol-gl/webgl/compat';
import { PixelData } from '../mol-util/image';
import { readTexture } from '../mol-gl/compute/util';
import { DrawPass } from './passes/draw';
import { PickPass } from './passes/pick';
import { Task } from '../mol-task';
import { ImagePass, ImageProps } from './passes/image';

export const Canvas3DParams = {
    cameraMode: PD.Select('perspective', [['perspective', 'Perspective'], ['orthographic', 'Orthographic']]),
    cameraFog: PD.Numeric(50, { min: 1, max: 100, step: 1 }),
    cameraResetDurationMs: PD.Numeric(250, { min: 0, max: 1000, step: 1 }, { description: 'The time it takes to reset the camera.' }),

    multiSample: PD.Group(MultiSampleParams),
    postprocessing: PD.Group(PostprocessingParams),
    renderer: PD.Group(RendererParams),
    trackball: PD.Group(TrackballControlsParams),
    debug: PD.Group(DebugHelperParams)
}
export const DefaultCanvas3DParams = PD.getDefaultValues(Canvas3DParams);
export type Canvas3DProps = PD.Values<typeof Canvas3DParams>

export { Canvas3D }

interface Canvas3D {
    readonly webgl: WebGLContext,

    add: (repr: Representation.Any) => Promise<void>
    remove: (repr: Representation.Any) => Promise<void>
    update: (repr?: Representation.Any, keepBoundingSphere?: boolean) => void
    clear: () => void

    // draw: (force?: boolean) => void
    requestDraw: (force?: boolean) => void
    animate: () => void
    identify: (x: number, y: number) => PickingId | undefined
    mark: (loci: Representation.Loci, action: MarkerAction) => void
    getLoci: (pickingId: PickingId) => Representation.Loci

    readonly didDraw: BehaviorSubject<now.Timestamp>
    readonly reprCount: BehaviorSubject<number>

    handleResize: () => void
    /** Focuses camera on scene's bounding sphere, centered and zoomed. */
    resetCamera: () => void
    readonly camera: Camera
    downloadScreenshot: () => void
    getPixelData: (variant: GraphicsRenderVariant) => PixelData
    setProps: (props: Partial<Canvas3DProps>) => void
    getImagePass: () => ImagePass

    /** Returns a copy of the current Canvas3D instance props */
    readonly props: Readonly<Canvas3DProps>
    readonly input: InputObserver
    readonly stats: RendererStats
    readonly interaction: Canvas3dInteractionHelper['events']

    dispose: () => void
}

const requestAnimationFrame = typeof window !== 'undefined' ? window.requestAnimationFrame : (f: (time: number) => void) => setImmediate(()=>f(Date.now()))
const DefaultRunTask = (task: Task<unknown>) => task.run()

namespace Canvas3D {
    export interface HoverEvent { current: Representation.Loci, buttons: ButtonsType, modifiers: ModifiersKeys }
    export interface ClickEvent { current: Representation.Loci, buttons: ButtonsType, modifiers: ModifiersKeys }

    export function fromCanvas(canvas: HTMLCanvasElement, props: Partial<Canvas3DProps> = {}, runTask = DefaultRunTask) {
        const gl = getGLContext(canvas, {
            alpha: true,
            antialias: true,
            depth: true,
            preserveDrawingBuffer: true,
            premultipliedAlpha: false,
        })
        if (gl === null) throw new Error('Could not create a WebGL rendering context')
        const input = InputObserver.fromElement(canvas)
        return Canvas3D.create(gl, input, props, runTask)
    }

    export function create(gl: GLRenderingContext, input: InputObserver, props: Partial<Canvas3DProps> = {}, runTask = DefaultRunTask): Canvas3D {
        const p = { ...DefaultCanvas3DParams, ...props }

        const reprRenderObjects = new Map<Representation.Any, Set<GraphicsRenderObject>>()
        const reprUpdatedSubscriptions = new Map<Representation.Any, Subscription>()
        const reprCount = new BehaviorSubject(0)

        const startTime = now()
        const didDraw = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp)

        const webgl = createContext(gl)

        let width = gl.drawingBufferWidth
        let height = gl.drawingBufferHeight

        const scene = Scene.create(webgl)

        const camera = new Camera({
            position: Vec3.create(0, 0, 100),
            mode: p.cameraMode,
            fog: p.cameraFog
        })

        const controls = TrackballControls.create(input, camera, p.trackball)
        const renderer = Renderer.create(webgl, p.renderer)
        const debugHelper = new BoundingSphereHelper(webgl, scene, p.debug);
        const interactionHelper = new Canvas3dInteractionHelper(identify, getLoci, input);

        const drawPass = new DrawPass(webgl, renderer, scene, camera, debugHelper)
        const pickPass = new PickPass(webgl, renderer, scene, camera, 0.5)
        const postprocessing = new PostprocessingPass(webgl, camera, drawPass, p.postprocessing)
        const multiSample = new MultiSamplePass(webgl, camera, drawPass, postprocessing, p.multiSample)

        let drawPending = false
        let cameraResetRequested = false

        function getLoci(pickingId: PickingId) {
            let loci: Loci = EmptyLoci
            let repr: Representation.Any = Representation.Empty
            reprRenderObjects.forEach((_, _repr) => {
                const _loci = _repr.getLoci(pickingId)
                if (!isEmptyLoci(_loci)) {
                    if (!isEmptyLoci(loci)) {
                        console.warn('found another loci, this should not happen')
                    }
                    loci = _loci
                    repr = _repr
                }
            })
            return { loci, repr }
        }

        function mark(reprLoci: Representation.Loci, action: MarkerAction) {
            const { repr, loci } = reprLoci
            let changed = false
            if (repr) {
                changed = repr.mark(loci, action)
            } else {
                reprRenderObjects.forEach((_, _repr) => { changed = _repr.mark(loci, action) || changed })
            }
            if (changed) {
                scene.update(void 0, true)
                const prevPickDirty = pickPass.pickDirty
                draw(true)
                pickPass.pickDirty = prevPickDirty // marking does not change picking buffers
            }
        }

        function render(variant: 'pick' | 'draw', force: boolean) {
            if (scene.isCommiting) return false

            let didRender = false
            controls.update(currentTime)
            Viewport.set(camera.viewport, 0, 0, width, height)
            const cameraChanged = camera.update()
            multiSample.update(force || cameraChanged, currentTime)

            if (force || cameraChanged || multiSample.enabled) {
                switch (variant) {
                    case 'pick':
                        pickPass.render()
                        break;
                    case 'draw':
                        renderer.setViewport(0, 0, width, height)
                        if (multiSample.enabled) {
                            multiSample.render(true)
                        } else {
                            drawPass.render(!postprocessing.enabled)
                            if (postprocessing.enabled) postprocessing.render(true)
                        }
                        pickPass.pickDirty = true
                        break;
                }
                didRender = true
            }

            return didRender;
        }

        let forceNextDraw = false;
        let currentTime = 0;

        function draw(force?: boolean) {
            if (render('draw', !!force || forceNextDraw)) {
                didDraw.next(now() - startTime as now.Timestamp)
            }
            forceNextDraw = false;
            drawPending = false
        }

        function requestDraw(force?: boolean) {
            if (drawPending) return
            drawPending = true
            forceNextDraw = !!force;
        }

        function animate() {
            currentTime = now();
            camera.transition.tick(currentTime);
            draw(false);
            if (!camera.transition.inTransition) interactionHelper.tick(currentTime);
            requestAnimationFrame(animate)
        }

        function identify(x: number, y: number): PickingId | undefined {
            return pickPass.identify(x, y)
        }

        function commit(renderObjects?: readonly GraphicsRenderObject[]) {
            scene.update(renderObjects, false)

            return runTask(scene.commit()).then(() => {
                if (cameraResetRequested && !scene.isCommiting) {
                    camera.focus(scene.boundingSphere.center, scene.boundingSphere.radius)
                    cameraResetRequested = false
                }
                if (debugHelper.isEnabled) debugHelper.update()
                requestDraw(true)
                reprCount.next(reprRenderObjects.size)
            })
        }

        function add(repr: Representation.Any) {
            const oldRO = reprRenderObjects.get(repr)
            const newRO = new Set<GraphicsRenderObject>()
            repr.renderObjects.forEach(o => newRO.add(o))

            if (oldRO) {
                if (!SetUtils.areEqual(newRO, oldRO)) {
                    for (const o of Array.from(newRO)) { if (!oldRO.has(o)) scene.add(o) }
                    for (const o of Array.from(oldRO)) { if (!newRO.has(o)) scene.remove(o) }
                }
            } else {
                repr.renderObjects.forEach(o => scene.add(o))
            }
            reprRenderObjects.set(repr, newRO)
            return commit(repr.renderObjects)
        }

        handleResize()

        return {
            webgl,

            add: (repr: Representation.Any) => {
                reprUpdatedSubscriptions.set(repr, repr.updated.subscribe(_ => {
                    if (!repr.state.syncManually) add(repr)
                }))
                return add(repr)
            },
            remove: (repr: Representation.Any) => {
                const updatedSubscription = reprUpdatedSubscriptions.get(repr)
                if (updatedSubscription) {
                    updatedSubscription.unsubscribe()
                }
                const renderObjects = reprRenderObjects.get(repr)
                if (renderObjects) {
                    renderObjects.forEach(o => scene.remove(o))
                    reprRenderObjects.delete(repr)
                    return commit()
                }
                return Promise.resolve()
            },
            update: (repr, keepSphere) => {
                if (repr) {
                    if (!reprRenderObjects.has(repr)) return;
                    scene.update(repr.renderObjects, !!keepSphere);
                } else {
                    scene.update(void 0, !!keepSphere)
                }
            },
            clear: () => {
                reprRenderObjects.clear()
                scene.clear()
                debugHelper.clear()
                requestDraw(true)
                reprCount.next(reprRenderObjects.size)
            },

            // draw,
            requestDraw,
            animate,
            identify,
            mark,
            getLoci,

            handleResize,
            resetCamera: () => {
                if (scene.isCommiting) {
                    cameraResetRequested = true
                } else {
                    camera.focus(scene.boundingSphere.center, scene.boundingSphere.radius, p.cameraResetDurationMs)
                    requestDraw(true);
                }
            },
            camera,
            downloadScreenshot: () => {
                // TODO
            },
            getPixelData: (variant: GraphicsRenderVariant) => {
                switch (variant) {
                    case 'color': return webgl.getDrawingBufferPixelData()
                    case 'pickObject': return pickPass.objectPickTarget.getPixelData()
                    case 'pickInstance': return pickPass.instancePickTarget.getPixelData()
                    case 'pickGroup': return pickPass.groupPickTarget.getPixelData()
                    case 'depth': return readTexture(webgl, drawPass.depthTexture) as PixelData
                }
            },
            didDraw,
            reprCount,
            setProps: (props: Partial<Canvas3DProps>) => {
                if (props.cameraMode !== undefined && props.cameraMode !== camera.state.mode) {
                    camera.setState({ mode: props.cameraMode })
                }
                if (props.cameraFog !== undefined && props.cameraFog !== camera.state.fog) {
                    camera.setState({ fog: props.cameraFog })
                }
                if (props.cameraResetDurationMs !== undefined) p.cameraResetDurationMs = props.cameraResetDurationMs

                if (props.postprocessing) postprocessing.setProps(props.postprocessing)
                if (props.multiSample) multiSample.setProps(props.multiSample)
                if (props.renderer) renderer.setProps(props.renderer)
                if (props.trackball) controls.setProps(props.trackball)
                if (props.debug) debugHelper.setProps(props.debug)
                requestDraw(true)
            },
            getImagePass: (props: Partial<ImageProps> = {}) => {
                return new ImagePass(webgl, renderer, scene, camera, debugHelper, props)
            },

            get props() {
                return {
                    cameraMode: camera.state.mode,
                    cameraFog: camera.state.fog,
                    cameraResetDurationMs: p.cameraResetDurationMs,

                    postprocessing: { ...postprocessing.props },
                    multiSample: { ...multiSample.props },
                    renderer: { ...renderer.props },
                    trackball: { ...controls.props },
                    debug: { ...debugHelper.props }
                }
            },
            get input() {
                return input
            },
            get stats() {
                return renderer.stats
            },
            get interaction() {
                return interactionHelper.events
            },
            dispose: () => {
                scene.clear()
                debugHelper.clear()
                input.dispose()
                controls.dispose()
                renderer.dispose()
                interactionHelper.dispose()
            }
        }

        function handleResize() {
            width = gl.drawingBufferWidth
            height = gl.drawingBufferHeight

            renderer.setViewport(0, 0, width, height)
            Viewport.set(camera.viewport, 0, 0, width, height)
            Viewport.set(controls.viewport, 0, 0, width, height)

            drawPass.setSize(width, height)
            pickPass.setSize(width, height)
            postprocessing.setSize(width, height)
            multiSample.setSize(width, height)

            requestDraw(true)
        }
    }
}