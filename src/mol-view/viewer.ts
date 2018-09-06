/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BehaviorSubject } from 'rxjs';

import { Vec3, Mat4, EPSILON } from 'mol-math/linear-algebra'
import InputObserver from 'mol-util/input/input-observer'
import * as SetUtils from 'mol-util/set'
import Renderer, { RendererStats } from 'mol-gl/renderer'
import { RenderObject } from 'mol-gl/render-object'

import TrackballControls from './controls/trackball'
import { Viewport } from './camera/util'
import { PerspectiveCamera } from './camera/perspective'
import { resizeCanvas } from './util';
import { createContext } from 'mol-gl/webgl/context';
import { Representation } from 'mol-geo/representation';
import { createRenderTarget } from 'mol-gl/webgl/render-target';
import Scene from 'mol-gl/scene';
import { RenderVariant } from 'mol-gl/webgl/render-item';
import { PickingId, decodeIdRGBA } from 'mol-geo/util/picking';
import { MarkerAction } from 'mol-geo/util/marker-data';
import { Loci, EmptyLoci, isEmptyLoci } from 'mol-model/loci';
import { Color } from 'mol-util/color';

interface Viewer {
    center: (p: Vec3) => void

    hide: (repr: Representation<any>) => void
    show: (repr: Representation<any>) => void

    add: (repr: Representation<any>) => void
    remove: (repr: Representation<any>) => void
    update: () => void
    clear: () => void

    draw: (force?: boolean) => void
    requestDraw: () => void
    animate: () => void
    pick: () => void
    identify: (x: number, y: number) => PickingId
    mark: (loci: Loci, action: MarkerAction) => void
    getLoci: (pickingId: PickingId) => Loci

    reprCount: BehaviorSubject<number>
    identified: BehaviorSubject<string>
    didDraw: BehaviorSubject<number>

    handleResize: () => void
    resetCamera: () => void
    downloadScreenshot: () => void
    getImageData: (variant: RenderVariant) => ImageData

    input: InputObserver
    stats: RendererStats
    dispose: () => void
}

function getWebGLContext(canvas: HTMLCanvasElement, contextAttributes?: WebGLContextAttributes) {
    function getContext(contextId: 'webgl' | 'experimental-webgl') {
        try {
            return canvas.getContext(contextId, contextAttributes)
        } catch (e) {
            return null
        }
    }
    return getContext('webgl') || getContext('experimental-webgl')
}

namespace Viewer {
    export function create(canvas: HTMLCanvasElement, container: Element): Viewer {
        const reprMap = new Map<Representation<any>, Set<RenderObject>>()
        const reprCount = new BehaviorSubject(0)
        const identified = new BehaviorSubject('')

        const startTime = performance.now()
        const didDraw = new BehaviorSubject(0)
        const input = InputObserver.create(canvas)

        const camera = PerspectiveCamera.create({
            near: 0.1,
            far: 10000,
            position: Vec3.create(0, 0, 50)
        })
        // camera.lookAt(Vec3.create(0, 0, 0))

        const gl = getWebGLContext(canvas, {
            alpha: false,
            antialias: true,
            depth: true,
            preserveDrawingBuffer: true
        })
        if (gl === null) {
            throw new Error('Could not create a WebGL rendering context')
        }
        const ctx = createContext(gl)

        const scene = Scene.create(ctx)
        // const controls = TrackballControls.create(input, scene, {})
        const controls = TrackballControls.create(input, camera, {})
        // const renderer = Renderer.create(ctx, camera, { clearColor: 0xFFFFFF })
        const renderer = Renderer.create(ctx, camera, { clearColor: Color(0x000000) })

        const pickScale = 1
        const pickWidth = Math.round(canvas.width * pickScale)
        const pickHeight = Math.round(canvas.height * pickScale)
        const objectPickTarget = createRenderTarget(ctx, pickWidth, pickHeight)
        const instancePickTarget = createRenderTarget(ctx, pickWidth, pickHeight)
        const groupPickTarget = createRenderTarget(ctx, pickWidth, pickHeight)

        let pickDirty = true
        let drawPending = false
        let lastRenderTime = -1
        const prevProjectionView = Mat4.zero()
        const prevSceneView = Mat4.zero()

        function getLoci(pickingId: PickingId) {
            let loci: Loci = EmptyLoci
            reprMap.forEach((_, repr) => {
                const _loci = repr.getLoci(pickingId)
                if (!isEmptyLoci(_loci)) {
                    if (!isEmptyLoci(loci)) console.warn('found another loci')
                    loci = _loci
                }
            })
            return loci
        }

        function mark(loci: Loci, action: MarkerAction) {
            // reprMap.forEach((roSet, repr) => repr.mark(loci, action))
            // scene.update()
            // requestDraw()
        }

        let nearPlaneDelta = 0
        function computeNearDistance() {
            const focusRadius = scene.boundingSphere.radius
            let dist = Vec3.distance(controls.target, camera.position)
            if (dist > focusRadius) return dist - focusRadius
            return 0
        }

        function render(variant: RenderVariant, force?: boolean) {
            // const p = scene.boundingSphere.center
            // console.log(p[0], p[1], p[2])
            // Vec3.set(controls.target, p[0], p[1], p[2])

            const focusRadius = scene.boundingSphere.radius
            const targetDistance = Vec3.distance(controls.target, camera.position)
            // console.log(targetDistance, controls.target, camera.position)
            let near = computeNearDistance() + nearPlaneDelta
            camera.near = Math.max(0.01, Math.min(near, targetDistance - 0.5))

            let fogNear = targetDistance - camera.near + 1 * focusRadius - nearPlaneDelta;
            let fogFar = targetDistance - camera.near + 2 * focusRadius - nearPlaneDelta;

            // console.log(fogNear, fogFar);
            camera.fogNear = Math.max(fogNear, 0.1);
            camera.fogFar = Math.max(fogFar, 0.2);

            // console.log(camera.fogNear, camera.fogFar, targetDistance)

            switch (variant) {
                case 'pickObject': objectPickTarget.bind(); break;
                case 'pickInstance': instancePickTarget.bind(); break;
                case 'pickGroup': groupPickTarget.bind(); break;
                case 'draw':
                    ctx.unbindFramebuffer();
                    renderer.setViewport(0, 0, canvas.width, canvas.height);
                    break;
            }
            let didRender = false
            controls.update()
            camera.update()
            if (force || !Mat4.areEqual(camera.projectionView, prevProjectionView, EPSILON.Value) || !Mat4.areEqual(scene.view, prevSceneView, EPSILON.Value)) {
                // console.log('foo', force, prevSceneView, scene.view)
                Mat4.copy(prevProjectionView, camera.projectionView)
                Mat4.copy(prevSceneView, scene.view)
                renderer.render(scene, variant)
                if (variant === 'draw') {
                    lastRenderTime = performance.now()
                    pickDirty = true
                    // pick()
                }
                didRender = true
            }
            return didRender
        }

        function draw(force?: boolean) {
            if (render('draw', force)) {
                didDraw.next(performance.now() - startTime)
            }
            drawPending = false
        }

        function requestDraw () {
            if (drawPending) return
            drawPending = true
            window.requestAnimationFrame(() => draw(true))
        }

        function animate () {
            draw(false)
            if (performance.now() - lastRenderTime > 500) {
                if (pickDirty) pick()
            }
            window.requestAnimationFrame(() => animate())
        }

        function pick() {
            console.log('pick')
            render('pickObject', pickDirty)
            render('pickInstance', pickDirty)
            render('pickGroup', pickDirty)

            pickDirty = false
        }

        function identify (x: number, y: number): PickingId {
            x *= ctx.pixelRatio
            y *= ctx.pixelRatio
            y = canvas.height - y // flip y

            const buffer = new Uint8Array(4)
            const xp = Math.round(x * pickScale)
            const yp = Math.round(y * pickScale)

            objectPickTarget.bind()
            ctx.readPixels(xp, yp, 1, 1, buffer)
            const objectId = decodeIdRGBA(buffer[0], buffer[1], buffer[2])

            instancePickTarget.bind()
            ctx.readPixels(xp, yp, 1, 1, buffer)
            const instanceId = decodeIdRGBA(buffer[0], buffer[1], buffer[2])

            groupPickTarget.bind()
            ctx.readPixels(xp, yp, 1, 1, buffer)
            const groupId = decodeIdRGBA(buffer[0], buffer[1], buffer[2])

            return { objectId, instanceId, groupId }
        }

        handleResize()

        return {
            center: (p: Vec3) => {
                Vec3.set(controls.target, p[0], p[1], p[2])
            },

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
                reprCount.next(reprMap.size)
                scene.update()
            },
            remove: (repr: Representation<any>) => {
                const renderObjectSet = reprMap.get(repr)
                if (renderObjectSet) {
                    renderObjectSet.forEach(o => scene.remove(o))
                    reprMap.delete(repr)
                    reprCount.next(reprMap.size)
                    scene.update()
                }
            },
            update: () => scene.update(),
            clear: () => {
                reprMap.clear()
                scene.clear()
            },

            draw,
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
            reprCount,
            identified,
            didDraw,

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

export default Viewer