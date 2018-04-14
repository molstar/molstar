/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Subject } from 'rxjs';

import { Vec2 } from 'mol-math/linear-algebra';

import MouseWheel from './mouse-wheel'
import TouchPinch from './touch-pinch'
import { eventOffset } from './event-offset'

export function getButtons(event: MouseEvent | Touch) {
    if (typeof event === 'object') {
        if ('buttons' in event) {
            return event.buttons
        } else if ('which' in event) {
            const b = (event as any).which  // 'any' to support older browsers
            if (b === 2) {
                return 4
            } else if (b === 3) {
                return 2
            } else if (b > 0) {
                return 1<<(b-1)
            }
        } else if ('button' in event) {
            const b = (event as any).button  // 'any' to support older browsers
            if (b === 1) {
                return 4
            } else if (b === 2) {
                return 2
            } else if (b >= 0) {
                return 1<<b
            }
        }
    }
    return 0
}

export const DefaultInputObserverProps = {
    parent: window as Window | Element,
    noScroll: true
}
export type InputObserverProps = Partial<typeof DefaultInputObserverProps>

export type MouseModifiers = {
    shift: boolean,
    alt: boolean,
    control: boolean,
    meta: boolean
}

interface InputObserver {
    noScroll: boolean
    isDragging: () => boolean
    isPinching: () => boolean

    drag: Subject<{ dx: number, dy: number, buttons: number, modifiers: MouseModifiers }>,
    wheel: Subject<{ dx: number, dy: number, dz: number, event: WheelEvent }>,
    pinch: Subject<number>,
    // click: Subject<{ x: number, y: number, buttons: number, modifiers: MouseModifiers }>,

    dispose: () => void
}

namespace InputObserver {
    export function create (element: Element, props: InputObserverProps = {}): InputObserver {
        const { parent, noScroll } = { ...DefaultInputObserverProps, ...props }

        const mouseStart = Vec2.zero()
        const tmp = Vec2.zero()
        const tmp2 = Vec2.zero()
        const modifiers: MouseModifiers = {
            shift: false,
            alt: false,
            control: false,
            meta: false
        }

        const touchPinch = TouchPinch.create(element)
        const mouseWheel = MouseWheel.create(element, noScroll)

        let dragging = false
        let disposed = false
        let buttons = 0

        const drag = new Subject<{ dx: number, dy: number, buttons: number, modifiers: MouseModifiers }>()
        const wheel = mouseWheel.wheel
        const pinch = new Subject<number>()

        attach()

        return {
            get noScroll () { return mouseWheel.noScroll },
            set noScroll (value: boolean) { mouseWheel.noScroll = value },
            isDragging: () => dragging,
            isPinching,

            drag,
            wheel,
            pinch,

            dispose
        }

        function attach () {
            element.addEventListener('mousedown', onInputDown as any, false)

            // for dragging to work outside canvas bounds,
            // mouse move/up events have to be added to parent, i.e. window
            parent.addEventListener('mousemove', onInputMove as any, false)
            parent.addEventListener('mouseup', onInputUp as any, false)

            // don't allow simulated mouse events
            element.addEventListener('touchstart', preventDefault as any, false)

            element.addEventListener('touchmove', onTouchMove as any, false)

            touchPinch.place.subscribe(onPinchPlace)
            touchPinch.lift.subscribe(onPinchLift)
            touchPinch.change.subscribe(onPinchChange)

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
            if (disposed) return
            disposed = true

            mouseWheel.dispose()
            touchPinch.dispose()

            element.removeEventListener('touchstart', preventDefault as any, false)
            element.removeEventListener('touchmove', onTouchMove as any, false)

            element.removeEventListener('mousedown', onInputDown as any, false)

            parent.removeEventListener('mousemove', onInputMove as any, false)
            parent.removeEventListener('mouseup', onInputUp as any, false)

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

        function preventDefault (ev: Event | Touch) {
            if ('preventDefault' in ev) ev.preventDefault()
        }

        function handleBlur () {
            if (buttons || modifiers.shift || modifiers.alt || modifiers.meta || modifiers.control) {
                buttons = 0
                modifiers.shift = modifiers.alt = modifiers.control = modifiers.meta = false
            }
        }

        function handleMods (event: MouseEvent | KeyboardEvent) {
            if ('altKey' in event) modifiers.alt = !!event.altKey
            if ('shiftKey' in event) modifiers.shift = !!event.shiftKey
            if ('ctrlKey' in event) modifiers.control = !!event.ctrlKey
            if ('metaKey' in event) modifiers.meta = !!event.metaKey
        }

        function onTouchMove (ev: TouchEvent) {
            if (!dragging || isPinching()) return

            // find currently active finger
            for (let i = 0; i < ev.changedTouches.length; i++) {
                const changed = ev.changedTouches[i]
                const idx = touchPinch.indexOfTouch(changed)
                if (idx !== -1) {
                    onInputMove(changed)
                    break
                }
            }
        }

        function onPinchPlace ({ newTouch, oldTouch }: { newTouch?: Touch, oldTouch?: Touch }) {
            dragging = !isPinching()
            if (dragging) {
                const firstFinger = oldTouch || newTouch
                if (firstFinger) onInputDown(firstFinger)
            }
        }

        function onPinchLift ({ removed, otherTouch }: { removed?: Touch, otherTouch?: Touch }) {
            // if either finger is down, consider it dragging
            const sum = touchPinch.fingers.reduce((sum, item) => sum + (item ? 1 : 0), 0)
            dragging = sum >= 1

            if (dragging && otherTouch) {
                eventOffset(mouseStart, otherTouch, element)
            }
        }

        function isPinching () {
            return touchPinch.pinching
        }

        function onPinchChange ({ currentDistance, lastDistance }: { currentDistance: number, lastDistance: number }) {
            pinch.next(currentDistance - lastDistance)
        }

        function onInputDown (ev: MouseEvent | Touch) {
            preventDefault(ev)
            eventOffset(mouseStart, ev, element)
            if (insideBounds(mouseStart)) {
                dragging = true
            }
        }

        function onInputUp () {
            dragging = false
        }

        function onInputMove (ev: MouseEvent | Touch) {
            buttons = getButtons(ev)
            const end = eventOffset(tmp, ev, element)
            if (pinch && isPinching()) {
                Vec2.copy(mouseStart, end)
                return
            }
            if (!dragging) return
            const rect = getClientSize(tmp2)
            const dx = (end[0] - mouseStart[0]) / rect[0]
            const dy = (end[1] - mouseStart[1]) / rect[1]
            drag.next({ dx, dy, buttons, modifiers })
            mouseStart[0] = end[0]
            mouseStart[1] = end[1]
        }

        function insideBounds (pos: Vec2) {
            if (element instanceof Window || element instanceof Document || element === document.body) {
                return true
            } else {
                const rect = element.getBoundingClientRect()
                return pos[0] >= 0 && pos[1] >= 0 && pos[0] < rect.width && pos[1] < rect.height
            }
        }

        function getClientSize (out: Vec2) {
            let source = element
            if (source instanceof Window || source instanceof Document || source === document.body) {
                source = document.documentElement
            }
            out[0] = source.clientWidth
            out[1] = source.clientHeight
            return out
        }
    }
}

export default InputObserver