/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/mikolalysenko/mouse-wheel,
 * copyright (c) 2015 Mikola Lysenko. MIT License
 */

import toPixels from './to-pixels'

export type MouseWheelCallback = (dx: number, dy: number, dz: number, event: MouseWheelEvent) => void

export default function mouseWheelListen(element: Element, callback: MouseWheelCallback, noScroll = false) {
    const lineHeight = toPixels('ex', element)
    const listener = function (event: MouseWheelEvent) {
        if (noScroll) {
            event.preventDefault()
        }
        const mode = event.deltaMode
        let dx = event.deltaX || 0
        let dy = event.deltaY || 0
        let dz = event.deltaZ || 0
        let scale = 1
        switch (mode) {
            case 1: scale = lineHeight; break
            case 2: scale = window.innerHeight; break
        }
        dx *= scale
        dy *= scale
        dz *= scale
        if (dx || dy || dz) {
            return callback(dx, dy, dz, event)
        }
    }
    element.addEventListener('wheel', listener)
    return listener
}