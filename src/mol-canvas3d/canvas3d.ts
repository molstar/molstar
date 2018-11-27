/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BehaviorSubject, Subscription } from 'rxjs';
import { now } from 'mol-util/now';

import { Vec3 } from 'mol-math/linear-algebra'
import InputObserver from 'mol-util/input/input-observer'
import * as SetUtils from 'mol-util/set'
import Renderer, { RendererStats } from 'mol-gl/renderer'
import { RenderObject } from 'mol-gl/render-object'

import TrackballControls from './controls/trackball'
import { Viewport } from './camera/util'
import { resizeCanvas } from './util';
import { createContext, getGLContext, WebGLContext } from 'mol-gl/webgl/context';
import { Representation } from 'mol-repr/representation';
import { createRenderTarget } from 'mol-gl/webgl/render-target';
import Scene from 'mol-gl/scene';
import { RenderVariant } from 'mol-gl/webgl/render-item';
import { PickingId, decodeIdRGB } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { Loci, EmptyLoci, isEmptyLoci } from 'mol-model/loci';
import { Color } from 'mol-util/color';
import { Camera } from './camera';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { BoundingSphereHelper } from './helper/bounding-sphere-helper';

export const Canvas3DParams = {
    // TODO: FPS cap?
    // maxFps: PD.Numeric(30),
    cameraMode: PD.Select('perspective', [['perspective', 'Perspective'], ['orthographic', 'Orthographic']]),
    backgroundColor: PD.Color(Color(0x000000)),
    clip: PD.Interval([1, 100], { min: 1, max: 100, step: 1 }),
    fog: PD.Interval([50, 100], { min: 1, max: 100, step: 1 }),
    pickingAlphaThreshold: PD.Numeric(0.5, { min: 0.0, max: 1.0, step: 0.01 }, { description: 'The minimum opacity value needed for an object to be pickable.' }),
    debug: PD.Group({
        showBoundingSpheres: PD.Boolean(false, { description: 'Show bounding spheres of render objects.' }),
    })
}
export type Canvas3DParams = PD.Values<typeof Canvas3DParams>

export { Canvas3D }

interface Canvas3D {
    readonly webgl: WebGLContext,

    hide: (repr: Representation.Any) => void
    show: (repr: Representation.Any) => void

    add: (repr: Representation.Any) => void
    remove: (repr: Representation.Any) => void
    update: () => void
    clear: () => void

    // draw: (force?: boolean) => void
    requestDraw: (force?: boolean) => void
    animate: () => void
    pick: () => void
    identify: (x: number, y: number) => Promise<PickingId | undefined>
    mark: (loci: Loci, action: MarkerAction, repr?: Representation.Any) => void
    getLoci: (pickingId: PickingId) => { loci: Loci, repr?: Representation.Any }

    readonly didDraw: BehaviorSubject<now.Timestamp>

    handleResize: () => void
    /** Focuses camera on scene's bounding sphere, centered and zoomed. */
    resetCamera: () => void
    readonly camera: Camera
    downloadScreenshot: () => void
    getImageData: (variant: RenderVariant) => ImageData
    setProps: (props: Partial<Canvas3DParams>) => void

    /** Returns a copy of the current Canvas3D instance props */
    readonly props: Canvas3DParams
    readonly input: InputObserver
    readonly stats: RendererStats
    dispose: () => void
}

namespace Canvas3D {
    export function create(canvas: HTMLCanvasElement, container: Element, props: Partial<Canvas3DParams> = {}): Canvas3D {
        const p = { ...PD.getDefaultValues(Canvas3DParams), ...props }

        const reprRenderObjects = new Map<Representation.Any, Set<RenderObject>>()
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

        const gl = getGLContext(canvas, {
            alpha: false,
            antialias: true,
            depth: true,
            preserveDrawingBuffer: true
        })
        if (gl === null) {
            throw new Error('Could not create a WebGL rendering context')
        }
        const webgl = createContext(gl)

        const scene = Scene.create(webgl)
        const controls = TrackballControls.create(input, camera, {})
        const renderer = Renderer.create(webgl, camera, { clearColor: p.backgroundColor })

        const pickScale = 1
        const pickWidth = Math.round(canvas.width * pickScale)
        const pickHeight = Math.round(canvas.height * pickScale)
        const objectPickTarget = createRenderTarget(webgl, pickWidth, pickHeight)
        const instancePickTarget = createRenderTarget(webgl, pickWidth, pickHeight)
        const groupPickTarget = createRenderTarget(webgl, pickWidth, pickHeight)

        let pickDirty = true
        let isPicking = false
        let isUpdating = false
        let drawPending = false
        let lastRenderTime = -1

        const boundingSphereHelper = new BoundingSphereHelper(scene, p.debug.showBoundingSpheres)

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

        function mark(loci: Loci, action: MarkerAction, repr?: Representation.Any) {
            let changed = false
            reprRenderObjects.forEach((_, _repr) => {
                if (!repr || repr === _repr) {
                    changed = _repr.mark(loci, action) || changed
                }
            })
            if (changed) {
                scene.update(true)
                draw(true)
                pickDirty = false // picking buffers should not have changed
            }
        }

        let currentNear = -1, currentFar = -1, currentFogNear = -1, currentFogFar = -1
        function setClipping() {
            const cDist = Vec3.distance(camera.state.position, camera.state.target)
            const bRadius = Math.max(10, scene.boundingSphere.radius)

            const nearFactor = (50 - p.clip[0]) / 50
            const farFactor = -(50 - p.clip[1]) / 50
            const near = cDist - (bRadius * nearFactor)
            const far = cDist + (bRadius * farFactor)

            const fogNearFactor = (50 - p.fog[0]) / 50
            const fogFarFactor = -(50 - p.fog[1]) / 50
            const fogNear = cDist - (bRadius * fogNearFactor)
            const fogFar = cDist + (bRadius * fogFarFactor)

            if (near !== currentNear || far !== currentFar || fogNear !== currentFogNear || fogFar !== currentFogFar) {
                camera.setState({ near, far, fogNear, fogFar })
                currentNear = near, currentFar = far, currentFogNear = fogNear, currentFogFar = fogFar
            }
        }

        function render(variant: RenderVariant, force: boolean) {
            if (isPicking || isUpdating) return false

            let didRender = false
            controls.update()
            // TODO: is this a good fix? Also, setClipping does not work if the user has manually set a clipping plane.
            if (!camera.transition.inTransition) setClipping();
            const cameraChanged = camera.updateMatrices();

            if (force || cameraChanged) {
                switch (variant) {
                    case 'pickObject': objectPickTarget.bind(); break;
                    case 'pickInstance': instancePickTarget.bind(); break;
                    case 'pickGroup': groupPickTarget.bind(); break;
                    case 'draw':
                        webgl.unbindFramebuffer();
                        renderer.setViewport(0, 0, canvas.width, canvas.height);
                        break;
                }

                renderer.render(scene, variant)
                if (variant === 'draw') {
                    lastRenderTime = now()
                    pickDirty = true
                }
                didRender = true
            }

            return didRender && cameraChanged;
        }

        let forceNextDraw = false;

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
            const t = now();
            camera.transition.tick(t);
            draw(false)
            if (t - lastRenderTime > 200) {
                if (pickDirty) pick()
            }
            window.requestAnimationFrame(animate)
        }

        function pick() {
            render('pickObject', pickDirty)
            render('pickInstance', pickDirty)
            render('pickGroup', pickDirty)
            webgl.gl.finish()

            pickDirty = false
        }

        async function identify(x: number, y: number): Promise<PickingId | undefined> {
            if (pickDirty || isPicking) return undefined

            isPicking = true

            x *= webgl.pixelRatio
            y *= webgl.pixelRatio
            y = canvas.height - y // flip y

            const buffer = new Uint8Array(4)
            const xp = Math.round(x * pickScale)
            const yp = Math.round(y * pickScale)

            objectPickTarget.bind()
            // TODO slow in Chrome, ok in FF; doesn't play well with gpu surface calc
            // await webgl.readPixelsAsync(xp, yp, 1, 1, buffer)
            webgl.readPixels(xp, yp, 1, 1, buffer)
            const objectId = decodeIdRGB(buffer[0], buffer[1], buffer[2])
            if (objectId === -1) return

            instancePickTarget.bind()
            // await webgl.readPixelsAsync(xp, yp, 1, 1, buffer)
            webgl.readPixels(xp, yp, 1, 1, buffer)
            const instanceId = decodeIdRGB(buffer[0], buffer[1], buffer[2])
            if (instanceId === -1) return

            groupPickTarget.bind()
            // await webgl.readPixelsAsync(xp, yp, 1, 1, buffer)
            webgl.readPixels(xp, yp, 1, 1, buffer)
            const groupId = decodeIdRGB(buffer[0], buffer[1], buffer[2])
            if (groupId === -1) return

            isPicking = false

            return { objectId, instanceId, groupId }
        }

        function add(repr: Representation.Any) {
            isUpdating = true
            const oldRO = reprRenderObjects.get(repr)
            const newRO = new Set<RenderObject>()
            repr.renderObjects.forEach(o => newRO.add(o))
            if (oldRO) {
                SetUtils.difference(newRO, oldRO).forEach(o => scene.add(o))
                SetUtils.difference(oldRO, newRO).forEach(o => scene.remove(o))
                scene.update()
            } else {
                repr.renderObjects.forEach(o => scene.add(o))
            }
            reprRenderObjects.set(repr, newRO)
            reprCount.next(reprRenderObjects.size)
            boundingSphereHelper.update()
            scene.update()
            isUpdating = false
            requestDraw(true)
        }

        handleResize()

        return {
            webgl,

            hide: (repr: Representation.Any) => {
                const renderObjectSet = reprRenderObjects.get(repr)
                if (renderObjectSet) renderObjectSet.forEach(o => o.state.visible = false)
            },
            show: (repr: Representation.Any) => {
                const renderObjectSet = reprRenderObjects.get(repr)
                if (renderObjectSet) renderObjectSet.forEach(o => o.state.visible = true)
            },

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
                    reprCount.next(reprRenderObjects.size)
                    boundingSphereHelper.update()
                    scene.update()
                    isUpdating = false
                    requestDraw(true)
                }
            },
            update: () => scene.update(),
            clear: () => {
                reprRenderObjects.clear()
                scene.clear()
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
            getImageData: (variant: RenderVariant) => {
                switch (variant) {
                    case 'draw': return renderer.getImageData()
                    case 'pickObject': return objectPickTarget.getImageData()
                    case 'pickInstance': return instancePickTarget.getImageData()
                    case 'pickGroup': return groupPickTarget.getImageData()
                }
            },
            didDraw,
            setProps: (props: Partial<Canvas3DParams>) => {
                if (props.cameraMode !== undefined && props.cameraMode !== camera.state.mode) {
                    camera.setState({ mode: props.cameraMode })
                }
                if (props.backgroundColor !== undefined && props.backgroundColor !== renderer.props.clearColor) {
                    renderer.setClearColor(props.backgroundColor)
                }

                if (props.clip !== undefined) p.clip = [props.clip[0], props.clip[1]]
                if (props.fog !== undefined) p.fog = [props.fog[0], props.fog[1]]

                if (props.pickingAlphaThreshold !== undefined && props.pickingAlphaThreshold !== renderer.props.pickingAlphaThreshold) {
                    renderer.setPickingAlphaThreshold(props.pickingAlphaThreshold)
                }
                if (props.debug && props.debug.showBoundingSpheres !== undefined) {
                    boundingSphereHelper.visible = props.debug.showBoundingSpheres
                }
                requestDraw(true)
            },

            get props() {
                return {
                    cameraMode: camera.state.mode,
                    backgroundColor: renderer.props.clearColor,
                    clip: p.clip,
                    fog: p.fog,
                    pickingAlphaThreshold: renderer.props.pickingAlphaThreshold,
                    debug: {
                        showBoundingSpheres: boundingSphereHelper.visible
                    }
                }
            },
            get input() {
                return input
            },
            get stats() {
                return renderer.stats
            },
            dispose: () => {
                scene.clear()
                input.dispose()
                controls.dispose()
                renderer.dispose()
                camera.dispose()
            }
        }

        function handleResize() {
            resizeCanvas(canvas, container)
            renderer.setViewport(0, 0, canvas.width, canvas.height)
            Viewport.set(camera.viewport, 0, 0, canvas.width, canvas.height)
            Viewport.set(controls.viewport, 0, 0, canvas.width, canvas.height)

            const pickWidth = Math.round(canvas.width * pickScale)
            const pickHeight = Math.round(canvas.height * pickScale)
            objectPickTarget.setSize(pickWidth, pickHeight)
            instancePickTarget.setSize(pickWidth, pickHeight)
            groupPickTarget.setSize(pickWidth, pickHeight)
        }
    }
}