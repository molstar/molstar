/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4, EPSILON } from 'mol-math/linear-algebra'
import InputObserver from 'mol-util/input/input-observer'
import * as SetUtils from 'mol-util/set'
import Renderer, { RendererStats } from 'mol-gl/renderer'
import { RenderObject } from 'mol-gl/scene'
import { StructureRepresentation } from 'mol-geo/representation/structure';

import TrackballControls from './controls/trackball'
import { Viewport } from './camera/util'
import { PerspectiveCamera } from './camera/perspective'
import { resizeCanvas } from './util';
import { createContext } from 'mol-gl/webgl/context';

interface Viewer {
    hide: (repr: StructureRepresentation) => void
    show: (repr: StructureRepresentation) => void

    add: (repr: StructureRepresentation) => void
    remove: (repr: StructureRepresentation) => void
    update: () => void
    clear: () => void

    draw: (force?: boolean) => void
    requestDraw: () => void
    animate: () => void

    handleResize: () => void

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
        const reprMap = new Map<StructureRepresentation, Set<RenderObject>>()

        const input = InputObserver.create(canvas)
        input.resize.subscribe(handleResize)

        const camera = PerspectiveCamera.create({
            near: 0.01,
            far: 10000,
            position: Vec3.create(0, 0, 50)
        })

        const controls = TrackballControls.create(input, camera, {

        })

        const gl = getWebGLContext(canvas)
        if (gl === null) {
            throw new Error('Could not create a WebGL rendering context')
        }
        const ctx = createContext(gl)

        const renderer = Renderer.create(ctx, camera)

        let drawPending = false
        const prevProjectionView = Mat4.zero()

        function draw (force?: boolean) {
            controls.update()
            camera.update()
            if (force || !Mat4.areEqual(camera.projectionView, prevProjectionView, EPSILON.Value)) {
                Mat4.copy(prevProjectionView, camera.projectionView)
                renderer.draw()
            }
            drawPending = false
        }

        function requestDraw () {
            if (drawPending) return
            drawPending = true
            window.requestAnimationFrame(() => draw(true))
        }

        function animate () {
            draw()
            window.requestAnimationFrame(() => animate())
        }

        handleResize()

        return {
            hide: (repr: StructureRepresentation) => {
                const renderObjectSet = reprMap.get(repr)
                if (renderObjectSet) renderObjectSet.forEach(o => o.visible = false)
            },
            show: (repr: StructureRepresentation) => {
                const renderObjectSet = reprMap.get(repr)
                if (renderObjectSet) renderObjectSet.forEach(o => o.visible = true)
            },

            add: (repr: StructureRepresentation) => {
                const oldRO = reprMap.get(repr)
                const newRO = new Set<RenderObject>()
                repr.renderObjects.forEach(o => newRO.add(o))
                if (oldRO) {
                    SetUtils.difference(newRO, oldRO).forEach(o => renderer.add(o))
                    SetUtils.difference(oldRO, newRO).forEach(o => renderer.remove(o))
                } else {
                    repr.renderObjects.forEach(o => renderer.add(o))
                }
                reprMap.set(repr, newRO)
            },
            remove: (repr: StructureRepresentation) => {
                const renderObjectSet = reprMap.get(repr)
                if (renderObjectSet) renderObjectSet.forEach(o => renderer.remove(o))
            },
            update: () => renderer.update(),
            clear: () => {
                reprMap.clear()
                renderer.clear()
            },

            draw,
            requestDraw,
            animate,

            handleResize,

            get stats() {
                return renderer.stats
            },
            dispose: () => {
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
        }
    }
}

export default Viewer