/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BehaviorSubject, Subscription } from 'rxjs';
import { now } from 'mol-util/now';

import { Vec3, Vec2 } from 'mol-math/linear-algebra'
import InputObserver, { ModifiersKeys, ButtonsType } from 'mol-util/input/input-observer'
import Renderer, { RendererStats, RendererParams } from 'mol-gl/renderer'
import { GraphicsRenderObject } from 'mol-gl/render-object'

import { TrackballControls, TrackballControlsParams } from './controls/trackball'
import { Viewport } from './camera/util'
import { resizeCanvas } from './util';
import { createContext, getGLContext, WebGLContext } from 'mol-gl/webgl/context';
import { Representation } from 'mol-repr/representation';
import { createRenderTarget } from 'mol-gl/webgl/render-target';
import Scene from 'mol-gl/scene';
import { GraphicsRenderVariant } from 'mol-gl/webgl/render-item';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { Loci, EmptyLoci, isEmptyLoci } from 'mol-model/loci';
import { Camera } from './camera';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { BoundingSphereHelper, DebugHelperParams } from './helper/bounding-sphere-helper';
import { decodeFloatRGB } from 'mol-util/float-packing';
import { SetUtils } from 'mol-util/set';
import { Canvas3dInteractionHelper } from './helper/interaction-events';
import { createTexture } from 'mol-gl/webgl/texture';
import { ValueCell } from 'mol-util';
import { getPostprocessingRenderable, PostprocessingParams, setPostprocessingProps } from './helper/postprocessing';
import { JitterVectors, getComposeRenderable } from './helper/multi-sample';

export const Canvas3DParams = {
    // TODO: FPS cap?
    // maxFps: PD.Numeric(30),
    cameraMode: PD.Select('perspective', [['perspective', 'Perspective'], ['orthographic', 'Orthographic']]),
    cameraClipDistance: PD.Numeric(0, { min: 0.0, max: 50.0, step: 0.1 }, { description: 'The distance between camera and scene at which to clip regardless of near clipping plane.' }),
    clip: PD.Interval([1, 100], { min: 1, max: 100, step: 1 }),
    fog: PD.Interval([50, 100], { min: 1, max: 100, step: 1 }),

    multiSample: PD.Select('off', [['off', 'Off'], ['on', 'On'], ['temporal', 'Temporal']]),
    sampleLevel: PD.Numeric(2, { min: 0, max: 5, step: 1 }),

    postprocessing: PD.Group(PostprocessingParams),
    renderer: PD.Group(RendererParams),
    trackball: PD.Group(TrackballControlsParams),
    debug: PD.Group(DebugHelperParams)
}
export type Canvas3DProps = PD.Values<typeof Canvas3DParams>

export { Canvas3D }

interface Canvas3D {
    readonly webgl: WebGLContext,

    add: (repr: Representation.Any) => void
    remove: (repr: Representation.Any) => void
    update: (repr?: Representation.Any, keepBoundingSphere?: boolean) => void
    clear: () => void

    // draw: (force?: boolean) => void
    requestDraw: (force?: boolean) => void
    animate: () => void
    pick: () => void
    identify: (x: number, y: number) => PickingId | undefined
    mark: (loci: Representation.Loci, action: MarkerAction) => void
    getLoci: (pickingId: PickingId) => Representation.Loci

    readonly didDraw: BehaviorSubject<now.Timestamp>

    handleResize: () => void
    /** Focuses camera on scene's bounding sphere, centered and zoomed. */
    resetCamera: () => void
    readonly camera: Camera
    downloadScreenshot: () => void
    getImageData: (variant: GraphicsRenderVariant) => ImageData
    setProps: (props: Partial<Canvas3DProps>) => void

    /** Returns a copy of the current Canvas3D instance props */
    readonly props: Readonly<Canvas3DProps>
    readonly input: InputObserver
    readonly stats: RendererStats
    readonly interaction: Canvas3dInteractionHelper['events']

    dispose: () => void
}

namespace Canvas3D {
    export interface HighlightEvent { current: Representation.Loci, prev: Representation.Loci, modifiers?: ModifiersKeys }
    export interface ClickEvent { current: Representation.Loci, buttons: ButtonsType, modifiers: ModifiersKeys }

    export function create(canvas: HTMLCanvasElement, container: Element, props: Partial<Canvas3DProps> = {}): Canvas3D {
        const p = { ...PD.getDefaultValues(Canvas3DParams), ...props }

        const reprRenderObjects = new Map<Representation.Any, Set<GraphicsRenderObject>>()
        const reprUpdatedSubscriptions = new Map<Representation.Any, Subscription>()
        const reprCount = new BehaviorSubject(0)

        const startTime = now()
        const didDraw = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp)
        const input = InputObserver.create(canvas)

        const camera = new Camera({
            near: 0.1,
            far: 10000,
            position: Vec3.create(0, 0, 10),
            mode: p.cameraMode
        })

        const _gl = getGLContext(canvas, {
            alpha: false,
            antialias: true,
            depth: true,
            preserveDrawingBuffer: true
        })
        if (_gl === null) {
            throw new Error('Could not create a WebGL rendering context')
        }
        const webgl = createContext(_gl)
        const { state, gl } = webgl

        const scene = Scene.create(webgl)
        const controls = TrackballControls.create(input, camera, p.trackball)
        const renderer = Renderer.create(webgl, camera, p.renderer)

        const drawTarget = createRenderTarget(webgl, canvas.width, canvas.height)
        const depthTexture = createTexture(webgl, 'image-depth', 'depth', 'ushort', 'nearest')
        depthTexture.define(canvas.width, canvas.height)
        depthTexture.attachFramebuffer(drawTarget.framebuffer, 'depth')

        const postprocessingTarget = createRenderTarget(webgl, canvas.width, canvas.height)
        const postprocessing = getPostprocessingRenderable(webgl, drawTarget.texture, depthTexture, p.postprocessing)

        const composeTarget = createRenderTarget(webgl, canvas.width, canvas.height)
        const holdTarget = createRenderTarget(webgl, canvas.width, canvas.height)
        const compose = getComposeRenderable(webgl, drawTarget.texture)

        const pickBaseScale = 0.25
        let pickScale = pickBaseScale / webgl.pixelRatio
        let pickWidth = Math.round(canvas.width * pickScale)
        let pickHeight = Math.round(canvas.height * pickScale)
        const objectPickTarget = createRenderTarget(webgl, pickWidth, pickHeight)
        const instancePickTarget = createRenderTarget(webgl, pickWidth, pickHeight)
        const groupPickTarget = createRenderTarget(webgl, pickWidth, pickHeight)

        let pickDirty = true
        let isIdentifying = false
        let isUpdating = false
        let drawPending = false

        let multiSampleIndex = -1
        let multiSample = false

        const debugHelper = new BoundingSphereHelper(webgl, scene, p.debug);
        const interactionHelper = new Canvas3dInteractionHelper(identify, getLoci, input);

        function getLoci(pickingId: PickingId) {
            let loci: Loci = EmptyLoci
            let repr: Representation.Any = Representation.Empty
            reprRenderObjects.forEach((_, _repr) => {
                const _loci = _repr.getLoci(pickingId)
                if (!isEmptyLoci(_loci)) {
                    if (!isEmptyLoci(loci)) console.warn('found another loci')
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
                const prevPickDirty = pickDirty
                draw(true)
                pickDirty = prevPickDirty // marking does not change picking buffers
            }
        }

        let currentNear = -1, currentFar = -1, currentFogNear = -1, currentFogFar = -1
        function setClipping() {
            const cDist = Vec3.distance(camera.state.position, camera.state.target)
            const bRadius = Math.max(10, scene.boundingSphere.radius)

            const nearFactor = (50 - p.clip[0]) / 50
            const farFactor = -(50 - p.clip[1]) / 50
            let near = cDist - (bRadius * nearFactor)
            let far = cDist + (bRadius * farFactor)

            const fogNearFactor = (50 - p.fog[0]) / 50
            const fogFarFactor = -(50 - p.fog[1]) / 50
            let fogNear = cDist - (bRadius * fogNearFactor)
            let fogFar = cDist + (bRadius * fogFarFactor)

            if (camera.state.mode === 'perspective') {
                near = Math.max(1, p.cameraClipDistance, near)
                far = Math.max(1, far)
                fogNear = Math.max(1, fogNear)
                fogFar = Math.max(1, fogFar)
            } else if (camera.state.mode === 'orthographic') {
                if (p.cameraClipDistance > 0) {
                    near = Math.max(p.cameraClipDistance, near)
                }
            }

            if (near !== currentNear || far !== currentFar || fogNear !== currentFogNear || fogFar !== currentFogFar) {
                camera.setState({ near, far, fogNear, fogFar })
                currentNear = near, currentFar = far, currentFogNear = fogNear, currentFogFar = fogFar
            }
        }

        function renderDraw() {
            renderer.setViewport(0, 0, canvas.width, canvas.height)
            renderer.render(scene, 'draw')
            if (debugHelper.isEnabled) {
                debugHelper.syncVisibility()
                renderer.render(debugHelper.scene, 'draw')
            }
        }

        function renderPostprocessing() {
            gl.viewport(0, 0, canvas.width, canvas.height)
            state.disable(gl.SCISSOR_TEST)
            state.disable(gl.BLEND)
            state.disable(gl.DEPTH_TEST)
            state.depthMask(false)
            postprocessing.render()
        }

        function renderTemporalMultiSample(postprocessingEnabled: boolean) {
            // based on the Multisample Anti-Aliasing Render Pass
            // contributed to three.js by bhouston / http://clara.io/
            //
            // This manual approach to MSAA re-renders the scene once for
            // each sample with camera jitter and accumulates the results.
            const offsetList = JitterVectors[ Math.max(0, Math.min(p.sampleLevel, 5)) ]

            if (multiSampleIndex === -1) return
            if (multiSampleIndex >= offsetList.length) {
                multiSampleIndex = -1
                return
            }

            const i = multiSampleIndex

            if (i === 0) {
                drawTarget.bind()
                renderDraw()
                if (postprocessingEnabled) {
                    postprocessingTarget.bind()
                    renderPostprocessing()
                }
                ValueCell.update(compose.values.uWeight, 1.0)
                ValueCell.update(compose.values.tColor, postprocessingEnabled ? postprocessingTarget.texture : drawTarget.texture)
                compose.update()

                holdTarget.bind()
                state.disable(gl.BLEND)
                compose.render()
            }

            const sampleWeight = 1.0 / offsetList.length

            camera.viewOffset.enabled = true
            ValueCell.update(compose.values.tColor, postprocessingEnabled ? postprocessingTarget.texture : drawTarget.texture)
            ValueCell.update(compose.values.uWeight, sampleWeight)
            compose.update()

            const { width, height } = drawTarget

            // render the scene multiple times, each slightly jitter offset
            // from the last and accumulate the results.
            const numSamplesPerFrame = Math.pow(2, p.sampleLevel)
            for (let i = 0; i < numSamplesPerFrame; ++i) {
                const offset = offsetList[multiSampleIndex]
                Camera.setViewOffset(camera.viewOffset, width, height, offset[0], offset[1], width, height)
                camera.updateMatrices()

                // render scene and optionally postprocess
                drawTarget.bind()
                renderDraw()
                if (postprocessingEnabled) {
                    postprocessingTarget.bind()
                    renderPostprocessing()
                }

                // compose rendered scene with compose target
                composeTarget.bind()
                state.enable(gl.BLEND)
                state.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD)
                state.blendFuncSeparate(gl.ONE, gl.ONE, gl.ONE, gl.ONE)
                state.disable(gl.DEPTH_TEST)
                state.disable(gl.SCISSOR_TEST)
                state.depthMask(false)
                if (multiSampleIndex === 0) {
                    webgl.state.clearColor(0, 0, 0, 0)
                    gl.clear(gl.COLOR_BUFFER_BIT)
                }
                compose.render()

                multiSampleIndex += 1
                if (multiSampleIndex >= offsetList.length ) break
            }

            const accumulationWeight = multiSampleIndex * sampleWeight
            if (accumulationWeight > 0) {
                ValueCell.update(compose.values.uWeight, 1.0)
                ValueCell.update(compose.values.tColor, composeTarget.texture)
                compose.update()
                webgl.unbindFramebuffer()
                gl.viewport(0, 0, canvas.width, canvas.height)
                state.disable(gl.BLEND)
                compose.render()
            }
            if (accumulationWeight < 1.0) {
                ValueCell.update(compose.values.uWeight, 1.0 - accumulationWeight)
                ValueCell.update(compose.values.tColor, holdTarget.texture)
                compose.update()
                webgl.unbindFramebuffer()
                gl.viewport(0, 0, canvas.width, canvas.height)
                if (accumulationWeight === 0) state.disable(gl.BLEND)
                else state.enable(gl.BLEND)
                compose.render()
            }

            camera.viewOffset.enabled = false
            camera.updateMatrices()
            if (multiSampleIndex >= offsetList.length) multiSampleIndex = -1
        }

        function renderMultiSample(postprocessingEnabled: boolean) {
            // based on the Multisample Anti-Aliasing Render Pass
            // contributed to three.js by bhouston / http://clara.io/
            //
            // This manual approach to MSAA re-renders the scene once for
            // each sample with camera jitter and accumulates the results.
            const offsetList = JitterVectors[ Math.max(0, Math.min(p.sampleLevel, 5)) ]

            const baseSampleWeight = 1.0 / offsetList.length
            const roundingRange = 1 / 32

            camera.viewOffset.enabled = true
            ValueCell.update(compose.values.tColor, postprocessingEnabled ? postprocessingTarget.texture : drawTarget.texture)
            compose.update()

            const { width, height } = drawTarget

            // render the scene multiple times, each slightly jitter offset
            // from the last and accumulate the results.
            for (let i = 0; i < offsetList.length; ++i) {
                const offset = offsetList[i]
                Camera.setViewOffset(camera.viewOffset, width, height, offset[0], offset[1], width, height)
                camera.updateMatrices()

                // the theory is that equal weights for each sample lead to an accumulation of rounding
                // errors. The following equation varies the sampleWeight per sample so that it is uniformly
                // distributed across a range of values whose rounding errors cancel each other out.
                const uniformCenteredDistribution = -0.5 + (i + 0.5) / offsetList.length
                const sampleWeight = baseSampleWeight + roundingRange * uniformCenteredDistribution
                ValueCell.update(compose.values.uWeight, sampleWeight)

                // render scene and optionally postprocess
                drawTarget.bind()
                renderDraw()
                if (postprocessingEnabled) {
                    postprocessingTarget.bind()
                    renderPostprocessing()
                }

                // compose rendered scene with compose target
                composeTarget.bind()
                gl.viewport(0, 0, canvas.width, canvas.height)
                state.enable(gl.BLEND)
                state.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD)
                state.blendFuncSeparate(gl.ONE, gl.ONE, gl.ONE, gl.ONE)
                state.disable(gl.DEPTH_TEST)
                state.disable(gl.SCISSOR_TEST)
                state.depthMask(false)
                if (i === 0) {
                    webgl.state.clearColor(0, 0, 0, 0)
                    gl.clear(gl.COLOR_BUFFER_BIT)
                }
                compose.render()
            }

            ValueCell.update(compose.values.uWeight, 1.0)
            ValueCell.update(compose.values.tColor, composeTarget.texture)
            compose.update()

            webgl.unbindFramebuffer()
            gl.viewport(0, 0, canvas.width, canvas.height)
            state.disable(gl.BLEND)
            compose.render()

            camera.viewOffset.enabled = false
            camera.updateMatrices()
        }

        function render(variant: 'pick' | 'draw', force: boolean) {
            if (isIdentifying || isUpdating) return false

            let didRender = false
            controls.update(currentTime);
            // TODO: is this a good fix? Also, setClipping does not work if the user has manually set a clipping plane.
            if (!camera.transition.inTransition) setClipping();
            const cameraChanged = camera.updateMatrices();
            if (force || cameraChanged) lastRenderTime = now()
            if (p.multiSample === 'temporal') {
                if (currentTime - lastRenderTime > 200) {
                    multiSample = multiSampleIndex !== -1
                } else {
                    multiSampleIndex = 0
                    multiSample = false
                }
            } else if (p.multiSample === 'on') {
                multiSample = true
            } else {
                multiSample = false
            }
            const postprocessingEnabled = p.postprocessing.occlusionEnable || p.postprocessing.outlineEnable

            if (force || cameraChanged || multiSample) {
                switch (variant) {
                    case 'pick':
                        renderer.setViewport(0, 0, pickWidth, pickHeight);
                        objectPickTarget.bind();
                        renderer.render(scene, 'pickObject');
                        instancePickTarget.bind();
                        renderer.render(scene, 'pickInstance');
                        groupPickTarget.bind();
                        renderer.render(scene, 'pickGroup');
                        break;
                    case 'draw':
                        renderer.setViewport(0, 0, canvas.width, canvas.height);
                        if (multiSample) {
                            if (p.multiSample === 'temporal') {
                                renderTemporalMultiSample(postprocessingEnabled)
                            } else {
                                renderMultiSample(postprocessingEnabled)
                            }
                        } else {
                            if (postprocessingEnabled) drawTarget.bind()
                            else webgl.unbindFramebuffer()
                            renderDraw()
                            if (postprocessingEnabled) {
                                webgl.unbindFramebuffer()
                                renderPostprocessing()
                            }
                        }
                        pickDirty = true
                        break;
                }
                didRender = true
            }

            return didRender && cameraChanged;
        }

        let forceNextDraw = false;
        let currentTime = 0;
        let lastRenderTime = 0;

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
            window.requestAnimationFrame(animate)
        }

        function pick() {
            if (pickDirty) {
                render('pick', true)
                pickDirty = false
            }
        }

        const readBuffer = new Uint8Array(4)
        function identify(x: number, y: number): PickingId | undefined {
            if (isIdentifying) return

            pick() // must be called before setting `isIdentifying = true`
            isIdentifying = true

            x *= webgl.pixelRatio
            y *= webgl.pixelRatio
            y = canvas.height - y // flip y

            const xp = Math.round(x * pickScale)
            const yp = Math.round(y * pickScale)

            objectPickTarget.bind()
            webgl.readPixels(xp, yp, 1, 1, readBuffer)
            const objectId = decodeFloatRGB(readBuffer[0], readBuffer[1], readBuffer[2])
            if (objectId === -1) { isIdentifying = false; return; }

            instancePickTarget.bind()
            webgl.readPixels(xp, yp, 1, 1, readBuffer)
            const instanceId = decodeFloatRGB(readBuffer[0], readBuffer[1], readBuffer[2])
            if (instanceId === -1) { isIdentifying = false; return; }

            groupPickTarget.bind()
            webgl.readPixels(xp, yp, 1, 1, readBuffer)
            const groupId = decodeFloatRGB(readBuffer[0], readBuffer[1], readBuffer[2])
            if (groupId === -1) { isIdentifying = false; return; }

            isIdentifying = false

            return { objectId, instanceId, groupId }
        }

        function add(repr: Representation.Any) {
            isUpdating = true
            const oldRO = reprRenderObjects.get(repr)
            const newRO = new Set<GraphicsRenderObject>()
            repr.renderObjects.forEach(o => newRO.add(o))

            if (oldRO) {
                if (!SetUtils.areEqual(newRO, oldRO)) {
                    for (const o of Array.from(newRO)) { if (!oldRO.has(o)) scene.add(o); }
                    for (const o of Array.from(oldRO)) { if (!newRO.has(o)) scene.remove(o) }
                }
            } else {
                repr.renderObjects.forEach(o => scene.add(o))
            }
            reprRenderObjects.set(repr, newRO)
            scene.update(repr.renderObjects, false)
            if (debugHelper.isEnabled) debugHelper.update()
            isUpdating = false
            requestDraw(true)
            reprCount.next(reprRenderObjects.size)
        }

        handleResize()

        return {
            webgl,

            add: (repr: Representation.Any) => {
                add(repr)
                reprUpdatedSubscriptions.set(repr, repr.updated.subscribe(_ => {
                    if (!repr.state.syncManually) add(repr)
                }))
            },
            remove: (repr: Representation.Any) => {
                const updatedSubscription = reprUpdatedSubscriptions.get(repr)
                if (updatedSubscription) {
                    updatedSubscription.unsubscribe()
                }
                const renderObjects = reprRenderObjects.get(repr)
                if (renderObjects) {
                    isUpdating = true
                    renderObjects.forEach(o => scene.remove(o))
                    reprRenderObjects.delete(repr)
                    scene.update(void 0, false)
                    if (debugHelper.isEnabled) debugHelper.update()
                    isUpdating = false
                    requestDraw(true)
                    reprCount.next(reprRenderObjects.size)
                }
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
            },

            // draw,
            requestDraw,
            animate,
            pick,
            identify,
            mark,
            getLoci,

            handleResize,
            resetCamera: () => {
                camera.focus(scene.boundingSphere.center, scene.boundingSphere.radius)
                requestDraw(true);
            },
            camera,
            downloadScreenshot: () => {
                // TODO
            },
            getImageData: (variant: GraphicsRenderVariant) => {
                switch (variant) {
                    case 'draw': return renderer.getImageData()
                    case 'pickObject': return objectPickTarget.getImageData()
                    case 'pickInstance': return instancePickTarget.getImageData()
                    case 'pickGroup': return groupPickTarget.getImageData()
                }
            },
            didDraw,
            setProps: (props: Partial<Canvas3DProps>) => {
                if (props.cameraMode !== undefined && props.cameraMode !== camera.state.mode) {
                    camera.setState({ mode: props.cameraMode })
                }
                if (props.cameraClipDistance !== undefined) p.cameraClipDistance = props.cameraClipDistance
                if (props.clip !== undefined) p.clip = [props.clip[0], props.clip[1]]
                if (props.fog !== undefined) p.fog = [props.fog[0], props.fog[1]]

                if (props.multiSample !== undefined) p.multiSample = props.multiSample
                if (props.sampleLevel !== undefined) p.sampleLevel = props.sampleLevel

                if (props.postprocessing) {
                    setPostprocessingProps(props.postprocessing, postprocessing, p.postprocessing, webgl)
                }

                if (props.renderer) renderer.setProps(props.renderer)
                if (props.trackball) controls.setProps(props.trackball)
                if (props.debug) debugHelper.setProps(props.debug)
                requestDraw(true)
            },

            get props() {
                return {
                    cameraMode: camera.state.mode,
                    cameraClipDistance: p.cameraClipDistance,
                    clip: p.clip,
                    fog: p.fog,

                    multiSample: p.multiSample,
                    sampleLevel: p.sampleLevel,

                    postprocessing: { ...p.postprocessing },
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
                camera.dispose()
                interactionHelper.dispose()
            }
        }

        function handleResize() {
            resizeCanvas(canvas, container)
            renderer.setViewport(0, 0, canvas.width, canvas.height)
            Viewport.set(camera.viewport, 0, 0, canvas.width, canvas.height)
            Viewport.set(controls.viewport, 0, 0, canvas.width, canvas.height)

            drawTarget.setSize(canvas.width, canvas.height)
            postprocessingTarget.setSize(canvas.width, canvas.height)
            composeTarget.setSize(canvas.width, canvas.height)
            holdTarget.setSize(canvas.width, canvas.height)
            depthTexture.define(canvas.width, canvas.height)
            ValueCell.update(postprocessing.values.uTexSize, Vec2.set(postprocessing.values.uTexSize.ref.value, canvas.width, canvas.height))
            ValueCell.update(compose.values.uTexSize, Vec2.set(compose.values.uTexSize.ref.value, canvas.width, canvas.height))

            pickScale = pickBaseScale / webgl.pixelRatio
            pickWidth = Math.round(canvas.width * pickScale)
            pickHeight = Math.round(canvas.height * pickScale)
            objectPickTarget.setSize(pickWidth, pickHeight)
            instancePickTarget.setSize(pickWidth, pickHeight)
            groupPickTarget.setSize(pickWidth, pickHeight)

            requestDraw(true)
        }
    }
}