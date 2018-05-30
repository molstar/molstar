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

    reprCount: BehaviorSubject<number>
    didDraw: BehaviorSubject<number>

    handleResize: () => void
    resetCamera: () => void
    downloadScreenshot: () => void
    getImageData: () => ImageData
    getPickImageData: () => ImageData

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

        const startTime = performance.now()
        const didDraw = new BehaviorSubject(0)

        const input = InputObserver.create(canvas)
        input.resize.subscribe(handleResize)

        const camera = PerspectiveCamera.create({
            near: 0.1,
            far: 10000,
            position: Vec3.create(0, 0, 50)
        })

        const controls = TrackballControls.create(input, camera, {

        })

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
        const renderer = Renderer.create(ctx, camera)

        const rtScale = 1 / 4
        const renderTarget = createRenderTarget(ctx, Math.round(canvas.width * rtScale), Math.round(canvas.height * rtScale))

        let drawPending = false
        const prevProjectionView = Mat4.zero()

        function render(pick: boolean, force?: boolean) {
            let didRender = false
            controls.update()
            camera.update()
            if (force || !Mat4.areEqual(camera.projectionView, prevProjectionView, EPSILON.Value)) {
                Mat4.copy(prevProjectionView, camera.projectionView)
                renderer.render(scene, pick)
                didRender = true
            }
            return didRender
        }

        function draw(force?: boolean) {
            ctx.unbindFramebuffer()
            const viewport = { x: 0, y: 0, width: canvas.width, height: canvas.height }
            renderer.setViewport(viewport)
            if (render(false, force)) {
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
            window.requestAnimationFrame(() => animate())
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
            },
            remove: (repr: Representation<any>) => {
                const renderObjectSet = reprMap.get(repr)
                if (renderObjectSet) {
                    renderObjectSet.forEach(o => scene.remove(o))
                    reprMap.delete(repr)
                    reprCount.next(reprMap.size)
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
            pick: () => {
                renderTarget.bind()
                render(true, true)
            },

            handleResize,
            resetCamera: () => {
                // TODO
            },
            downloadScreenshot: () => {
                // TODO
            },
            getImageData: () => {
                return renderer.getImageData()
            },
            getPickImageData: () => {
                return renderTarget.getImageData()
            },
            reprCount,
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
            const viewport = { x: 0, y: 0, width: canvas.width, height: canvas.height }
            renderer.setViewport(viewport)
            Viewport.copy(camera.viewport, viewport)
            Viewport.copy(controls.viewport, viewport)
            renderTarget.setSize(Math.round(canvas.width * rtScale), Math.round(canvas.height * rtScale))
        }
    }
}

export default Viewer