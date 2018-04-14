/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/mikolalysenko/mouse-change,
 * copyright (c) 2015 Mikola Lysenko. MIT License
 */

import { Subject } from 'rxjs';

import * as mouse from './mouse-event'



interface MouseChange {
    change: Subject<number>,
    dispose: () => void
}

namespace MouseChange {
    export type Modifiers = {
        shift: boolean,
        alt: boolean,
        control: boolean,
        meta: boolean
    }
    export type Info = {
        buttons: number,
        x: number,
        y: number,
        modifiers: Modifiers
    }

    export function create(element: Element) {
        let buttonState = 0
        let x = 0
        let y = 0
        const mods: Modifiers = {
            shift: false,
            alt: false,
            control: false,
            meta: false
        }
        let attached = false

        const change = new Subject<Info>()

        // Attach listeners
        attachListeners()

        return {
            change,
            dispose
        }

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

                change.next({ buttons: buttonState, x, y, modifiers: mods })
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
                change.next({ buttons: 0, x: 0, y: 0, modifiers: mods })
            }
        }

        function handleMods (event: MouseEvent | KeyboardEvent) {
            if (updateMods(event)) {
                change.next({ buttons: buttonState, x, y, modifiers: mods })
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

        function dispose () {
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
    }
}

export default MouseChange