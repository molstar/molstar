/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** resize canvas to container element */
export function resizeCanvas (canvas: HTMLCanvasElement, container: Element) {
    let w = window.innerWidth
    let h = window.innerHeight
    if (container !== document.body) {
        let bounds = container.getBoundingClientRect()
        w = bounds.right - bounds.left
        h = bounds.bottom - bounds.top
    }
    canvas.width = window.devicePixelRatio * w
    canvas.height = window.devicePixelRatio * h
    Object.assign(canvas.style, { width: `${w}px`, height: `${h}px` })
}