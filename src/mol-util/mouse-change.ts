/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/mikolalysenko/mouse-change,
 * copyright (c) 2015 Mikola Lysenko. MIT License
 */

import * as mouse from './mouse-event'

export type MouseModifiers = {
    shift: boolean,
    alt: boolean,
    control: boolean,
    meta: boolean
}
export type MouseChangeCallback = (buttonState: number, x: number, y: number, mods: MouseModifiers) => void

export default function mouseListen (element: Element, callback: MouseChangeCallback) {
    let buttonState = 0
    let x = 0
    let y = 0
    const mods: MouseModifiers = {
        shift: false,
        alt: false,
        control: false,
        meta: false
    }
    let attached = false

    function updateMods (event: MouseEvent | KeyboardEvent) {
        let changed = false
        if ('altKey' in event) {
            changed = changed || event.altKey !== mods.alt
            mods.alt = !!event.altKey
        }
        if ('shiftKey' in event) {
            changed = changed || event.shiftKey !== mods.shift
            mods.shift = !!event.shiftKey
        }
        if ('ctrlKey' in event) {
            changed = changed || event.ctrlKey !== mods.control
            mods.control = !!event.ctrlKey
        }
        if ('metaKey' in event) {
            changed = changed || event.metaKey !== mods.meta
            mods.meta = !!event.metaKey
        }
        return changed
    }

    function handleEvent (nextButtons: number, event: MouseEvent) {
        const nextX = mouse.x(event)
        const nextY = mouse.y(event)
        if ('buttons' in event) {
            nextButtons = event.buttons | 0
        }
        if (nextButtons !== buttonState || nextX !== x || nextY !== y || updateMods(event) ) {
            buttonState = nextButtons | 0
            x = nextX || 0
            y = nextY || 0
            callback && callback(buttonState, x, y, mods)
        }
    }

    function clearState (event: MouseEvent) {
        handleEvent(0, event)
    }

    function handleBlur () {
        if (buttonState || x || y || mods.shift || mods.alt || mods.meta || mods.control) {
            x = y = 0
            buttonState = 0
            mods.shift = mods.alt = mods.control = mods.meta = false
            callback && callback(0, 0, 0, mods)
        }
    }

    function handleMods (event: MouseEvent | KeyboardEvent) {
        if (updateMods(event)) {
            callback && callback(buttonState, x, y, mods)
        }
    }

    function handleMouseMove (event: MouseEvent) {
        if (mouse.buttons(event) === 0) {
            handleEvent(0, event)
        } else {
            handleEvent(buttonState, event)
        }
    }

    function handleMouseDown (event: MouseEvent) {
        handleEvent(buttonState | mouse.buttons(event), event)
    }

    function handleMouseUp (event: MouseEvent) {
        handleEvent(buttonState & ~mouse.buttons(event), event)
    }

    function attachListeners () {
        if (attached) return
        attached = true

        element.addEventListener('mousemove', handleMouseMove as EventListener)
        element.addEventListener('mousedown', handleMouseDown as EventListener)
        element.addEventListener('mouseup', handleMouseUp as EventListener)

        element.addEventListener('mouseleave', clearState as EventListener)
        element.addEventListener('mouseenter', clearState as EventListener)
        element.addEventListener('mouseout', clearState as EventListener)
        element.addEventListener('mouseover', clearState as EventListener)

        element.addEventListener('blur', handleBlur)
        element.addEventListener('keyup', handleMods as EventListener)
        element.addEventListener('keydown', handleMods as EventListener)
        element.addEventListener('keypress', handleMods as EventListener)

        if (!(element instanceof Window)) {
            window.addEventListener('blur', handleBlur)
            window.addEventListener('keyup', handleMods)
            window.addEventListener('keydown', handleMods)
            window.addEventListener('keypress', handleMods)
        }
    }

    function detachListeners () {
        if (!attached) return
        attached = false

        element.removeEventListener('mousemove', handleMouseMove as EventListener)
        element.removeEventListener('mousedown', handleMouseDown as EventListener)
        element.removeEventListener('mouseup', handleMouseUp as EventListener)

        element.removeEventListener('mouseleave', clearState as EventListener)
        element.removeEventListener('mouseenter', clearState as EventListener)
        element.removeEventListener('mouseout', clearState as EventListener)
        element.removeEventListener('mouseover', clearState as EventListener)

        element.removeEventListener('blur', handleBlur)
        element.removeEventListener('keyup', handleMods as EventListener)
        element.removeEventListener('keydown', handleMods as EventListener)
        element.removeEventListener('keypress', handleMods as EventListener)

        if (!(element instanceof Window)) {
            window.removeEventListener('blur', handleBlur)
            window.removeEventListener('keyup', handleMods)
            window.removeEventListener('keydown', handleMods)
            window.removeEventListener('keypress', handleMods)
        }
    }

    // Attach listeners
    attachListeners()

    const result = {
        element: element
    }

    Object.defineProperties(result, {
        enabled: {
            get: function () { return attached },
            set: function (f) {
                if (f) {
                    attachListeners()
                } else {
                    detachListeners()
                }
            },
            enumerable: true
        },
        buttons: {
            get: function () { return buttonState },
            enumerable: true
        },
        x: {
            get: function () { return x },
            enumerable: true
        },
        y: {
            get: function () { return y },
            enumerable: true
        },
        mods: {
            get: function () { return mods },
            enumerable: true
        }
    })

    return result
}