/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/regl-project/regl-camera,
 * copyright (c) 2016 Mikola Lysenko. MIT License
 */

const isBrowser = typeof window !== 'undefined'

import REGL = require('regl');

import mouseChange, { MouseModifiers } from 'mol-util/mouse-change'
import mouseWheel from 'mol-util/mouse-wheel'
import { defaults } from 'mol-util'
import { Mat4, Vec3 } from 'mol-math/linear-algebra/3d'
import { clamp, damp } from 'mol-math/interpolate'

export interface CameraUniforms {
    projection: Mat4,
}

export interface CameraState {
    center: Vec3,
    theta: number,
    phi: number,
    distance: number,
    eye: Vec3,
    up: Vec3,
    fovy: number,
    near: number,
    far: number,
    noScroll: boolean,
    flipY: boolean,
    dtheta: number,
    dphi: number,
    rotationSpeed: number,
    zoomSpeed: number,
    renderOnDirty: boolean,
    damping: number,
    minDistance: number,
    maxDistance: number,
}

export interface Camera {
    update: (props: any, block: any) => void,
    setState: (newState: CameraState) => void,
    isDirty: () => boolean
}

export namespace Camera {
    export function create (regl: REGL.Regl, element: HTMLElement, props: Partial<CameraState> = {}): Camera {
        const state: CameraState = {
            center: defaults(props.center, Vec3.zero()),
            theta: defaults(props.theta, 0),
            phi: defaults(props.phi, 0),
            distance: Math.log(defaults(props.distance, 10.0)),
            eye: Vec3.zero(),
            up: defaults(props.up, Vec3.create(0, 1, 0)),
            fovy: defaults(props.fovy, Math.PI / 4.0),
            near: defaults(props.near, 0.01),
            far: defaults(props.far, 1000.0),
            noScroll: defaults(props.noScroll, false),
            flipY: defaults(props.flipY, false),
            dtheta: 0,
            dphi: 0,
            rotationSpeed: defaults(props.rotationSpeed, 1),
            zoomSpeed: defaults(props.zoomSpeed, 1),
            renderOnDirty: defaults(props.renderOnDirty, false),
            damping: defaults(props.damping, 0.9),
            minDistance: Math.log(defaults(props.minDistance, 0.1)),
            maxDistance: Math.log(defaults(props.maxDistance, 1000))
        }

        const view = Mat4.identity()
        const projection = Mat4.identity()

        const right = Vec3.create(1, 0, 0)
        const front = Vec3.create(0, 0, 1)

        let dirty = false
        let ddistance = 0

        let prevX = 0
        let prevY = 0

        if (isBrowser) {
            const source = element || regl._gl.canvas

            const getWidth = function () {
                return element ? element.offsetWidth : window.innerWidth
            }

            const getHeight = function () {
                return element ? element.offsetHeight : window.innerHeight
            }

            mouseChange(source, function (buttons: number, x: number, y: number, mods: MouseModifiers) {
                if (buttons & 1) {
                    const dx = (x - prevX) / getWidth()
                    const dy = (y - prevY) / getHeight()

                    state.dtheta += state.rotationSpeed * 4.0 * dx
                    state.dphi += state.rotationSpeed * 4.0 * dy
                    dirty = true;
                }
                prevX = x
                prevY = y
            })

            mouseWheel(source, function (dx: number, dy: number) {
                ddistance += dy / getHeight() * state.zoomSpeed
                dirty = true;
            }, state.noScroll)
        }

        function dampAndMarkDirty (x: number) {
            const xd = damp(x, state.damping)
            if (Math.abs(xd) < 0.1) return 0
            dirty = true;
            return xd
        }

        function setState (newState: Partial<CameraState> = {}) {
            Object.assign(state, newState)

            const { center, eye, up, dtheta, dphi } = state

            state.theta += dtheta
            state.phi = clamp(state.phi + dphi, -Math.PI / 2.0, Math.PI / 2.0)
            state.distance = clamp(state.distance + ddistance, state.minDistance, state.maxDistance)

            state.dtheta = dampAndMarkDirty(dtheta)
            state.dphi = dampAndMarkDirty(dphi)
            ddistance = dampAndMarkDirty(ddistance)

            const theta = state.theta
            const phi = state.phi
            const r = Math.exp(state.distance)

            const vf = r * Math.sin(theta) * Math.cos(phi)
            const vr = r * Math.cos(theta) * Math.cos(phi)
            const vu = r * Math.sin(phi)

            for (let i = 0; i < 3; ++i) {
                eye[i] = center[i] + vf * front[i] + vr * right[i] + vu * up[i]
            }

            Mat4.lookAt(view, eye, center, up)
        }

        const injectContext = regl({
            context: {
                view: () => view,
                dirty: () => dirty,
                projection: (context: REGL.DefaultContext) => {
                    Mat4.perspective(
                        projection,
                        state.fovy,
                        context.viewportWidth / context.viewportHeight,
                        state.near,
                        state.far
                    )
                    if (state.flipY) { projection[5] *= -1 }
                    return projection
                }
            },
            uniforms: {  // TODO
                view: regl.context('view' as any),
                projection: regl.context('projection' as any)
            }
        })

        function update (props: any, block: any) {
            setState()
            injectContext(props, block)
            if (dirty) {
                console.log(view)
            }
            dirty = false
        }

        return {
            update,
            setState,
            isDirty: () => dirty
        }
    }
}
