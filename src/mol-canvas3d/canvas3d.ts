/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BehaviorSubject } from 'rxjs';
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
import { Representation } from 'mol-repr';
import { createRenderTarget } from 'mol-gl/webgl/render-target';
import Scene from 'mol-gl/scene';
import { RenderVariant } from 'mol-gl/webgl/render-item';
import { PickingId, decodeIdRGB } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { Loci, EmptyLoci, isEmptyLoci } from 'mol-model/loci';
import { Color } from 'mol-util/color';
import { Camera } from './camera';

export const DefaultCanvas3DProps = {
    // TODO: FPS cap?
    // maxFps: 30,
    cameraPosition: Vec3.create(0, 0, 50),
    cameraMode: 'perspective' as Camera.Mode,
    backgroundColor: Color(0x000000),
}
export type Canvas3DProps = typeof DefaultCanvas3DProps

export { Canvas3D }

interface Canvas3D {
    readonly webgl: WebGLContext,

    hide: (repr: Representation<any>) => void
    show: (repr: Representation<any>) => void

    add: (repr: Representation<any>) => void
    remove: (repr: Representation<any>) => void
    update: () => void
    clear: () => void

    // draw: (force?: boolean) => void
    requestDraw: (force?: boolean) => void
    animate: () => void
    pick: () => void
    identify: (x: number, y: number) => Promise<PickingId | undefined>
    mark: (loci: Loci, action: MarkerAction) => void
    getLoci: (pickingId: PickingId) => { loci: Loci, repr?: Representation<any> }

    readonly didDraw: BehaviorSubject<now.Timestamp>

    handleResize: () => void
    resetCamera: () => void
    readonly camera: Camera
    downloadScreenshot: () => void
    getImageData: (variant: RenderVariant) => ImageData
    setProps: (props: Partial<Canvas3DProps>) => void

    /** Returns a copy of the current Canvas3D instance props */
    readonly props: Canvas3DProps
    readonly input: InputObserver
    readonly stats: RendererStats
    dispose: () => void
}

namespace Canvas3D {
    export function create(canvas: HTMLCanvasElement, container: Element, props: Partial<Canvas3DProps> = {}): Canvas3D {
        const p = { ...props, ...DefaultCanvas3DProps }

        const reprMap = new Map<Representation<any>, Set<RenderObject>>()

        const startTime = now()
        const didDraw = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp)
        const input = InputObserver.create(canvas)

        const camera = new Camera({
            near: 0.1,
            far: 10000,
            position: Vec3.clone(p.cameraPosition),
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
        let drawPending = false
        let lastRenderTime = -1

        function getLoci(pickingId: PickingId) {
            let loci: Loci = EmptyLoci
            let repr: Representation.Any = Representation.Empty
            reprMap.forEach((_, _repr) => {
                const _loci = _repr.getLoci(pickingId)
                if (!isEmptyLoci(_loci)) {
                    if (!isEmptyLoci(loci)) console.warn('found another loci')
                    loci = _loci
                    repr = _repr
                }
            })
            return { loci, repr }
        }

        function mark(loci: Loci, action: MarkerAction) {
            let changed = false
            reprMap.forEach((roSet, repr) => {
                changed = repr.mark(loci, action) || changed
            })
            if (changed) {
                // console.log('changed')
                scene.update()
                draw(true)
                pickDirty = false // picking buffers should not have changed
            }
        }

        // let nearPlaneDelta = 0
        // function computeNearDistance() {
        //     const focusRadius = scene.boundingSphere.radius
        //     let dist = Vec3.distance(controls.target, camera.position)
        //     if (dist > focusRadius) return dist - focusRadius
        //     return 0
        // }

        function render(variant: RenderVariant, force: boolean) {
            if (isPicking) return false
            // const p = scene.boundingSphere.center
            // console.log(p[0], p[1], p[2])
            // Vec3.set(controls.target, p[0], p[1], p[2])

            // TODO update near/far
            // const focusRadius = scene.boundingSphere.radius
            // const targetDistance = Vec3.distance(controls.target, camera.position)
            // console.log(targetDistance, controls.target, camera.position)
            // let near = computeNearDistance() + nearPlaneDelta
            // camera.near = Math.max(0.01, Math.min(near, targetDistance - 0.5))

            // let fogNear = targetDistance - camera.near + 1 * focusRadius - nearPlaneDelta;
            // let fogFar = targetDistance - camera.near + 2 * focusRadius - nearPlaneDelta;

            // // console.log(fogNear, fogFar);
            // camera.fogNear = Math.max(fogNear, 0.1);
            // camera.fogFar = Math.max(fogFar, 0.2);

            // console.log(camera.fogNear, camera.fogFar, targetDistance)

            let didRender = false
            controls.update()
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
            // The animation frame is being requested by animate already.
            // window.requestAnimationFrame(() => draw(force))
        }

        function animate() {
            draw(false)
            if (now() - lastRenderTime > 200) {
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
            if (pickDirty) return undefined

            isPicking = true

            x *= webgl.pixelRatio
            y *= webgl.pixelRatio
            y = canvas.height - y // flip y

            const buffer = new Uint8Array(4)
            const xp = Math.round(x * pickScale)
            const yp = Math.round(y * pickScale)

            objectPickTarget.bind()
            await webgl.readPixelsAsync(xp, yp, 1, 1, buffer)
            const objectId = decodeIdRGB(buffer[0], buffer[1], buffer[2])

            instancePickTarget.bind()
            await webgl.readPixels(xp, yp, 1, 1, buffer)
            const instanceId = decodeIdRGB(buffer[0], buffer[1], buffer[2])

            groupPickTarget.bind()
            await webgl.readPixels(xp, yp, 1, 1, buffer)
            const groupId = decodeIdRGB(buffer[0], buffer[1], buffer[2])

            isPicking = false

            // TODO
            if (objectId === -1 || instanceId === -1 || groupId === -1) {
                return { objectId: -1, instanceId: -1, groupId: -1 }
            } else {
                return { objectId, instanceId, groupId }
            }
        }

        handleResize()

        return {
            webgl,

            hide: (repr: Representation<any>) => {
                const renderObjectSet = reprMap.get(repr)
                if (renderObjectSet) renderObjectSet.forEach(o => o.state.visible = false)
            },
            show: (repr: Representation<any>) => {
                const renderObjectSet = reprMap.get(repr)
                if (renderObjectSet) renderObjectSet.forEach(o => o.state.visible = true)
            },

            add: (repr: Representation<any>) => {
                const oldRO = reprMap.get(repr)
                const newRO = new Set<RenderObject>()
                repr.renderObjects.forEach(o => newRO.add(o))
                if (oldRO) {
                    SetUtils.difference(newRO, oldRO).forEach(o => scene.add(o))
                    SetUtils.difference(oldRO, newRO).forEach(o => scene.remove(o))
                    scene.update()
                } else {
                    repr.renderObjects.forEach(o => scene.add(o))
                }
                reprMap.set(repr, newRO)
                scene.update()
            },
            remove: (repr: Representation<any>) => {
                const renderObjectSet = reprMap.get(repr)
                if (renderObjectSet) {
                    renderObjectSet.forEach(o => scene.remove(o))
                    reprMap.delete(repr)
                    scene.update()
                }
            },
            update: () => scene.update(),
            clear: () => {
                reprMap.clear()
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
                // TODO
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
            setProps: (props: Partial<Canvas3DProps>) => {
                if (props.cameraMode !== undefined && props.cameraMode !== camera.state.mode) {
                    camera.setState({ mode: props.cameraMode })
                }
                if (props.backgroundColor !== undefined && props.backgroundColor !== renderer.props.clearColor) {
                    renderer.setClearColor(props.backgroundColor)
                }
                requestDraw(true)
            },

            get props() {
                return {
                    cameraPosition: Vec3.clone(camera.position),
                    cameraMode: camera.state.mode,
                    backgroundColor: renderer.props.clearColor
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