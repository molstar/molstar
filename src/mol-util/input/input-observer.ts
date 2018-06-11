/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Subject } from 'rxjs';

import { Vec2 } from 'mol-math/linear-algebra';

import toPixels from '../to-pixels'

function getButtons(event: MouseEvent | Touch) {
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
    noScroll: true,
    noContextMenu: true
}
export type InputObserverProps = Partial<typeof DefaultInputObserverProps>

export type ModifiersKeys = {
    shift: boolean,
    alt: boolean,
    control: boolean,
    meta: boolean
}

export const enum ButtonsFlag {
    /** No button or un-initialized */
    None = 0x0,
    /** Primary button (usually left) */
    Primary = 0x1,
    /** Secondary button (usually right) */
    Secondary = 0x2,
    /** Auxilary button (usually middle or mouse wheel button)  */
    Auxilary = 0x4,
    /** 4th button (typically the "Browser Back" button) */
    Forth = 0x8,
    /** 5th button (typically the "Browser Forward" button) */
    Five = 0x10,
}

type BaseInput = {
    buttons: number
    modifiers: ModifiersKeys
}

export type DragInput = {
    x: number,
    y: number,
    dx: number,
    dy: number,
    pageX: number,
    pageY: number,
    isStart: boolean
} & BaseInput

export type WheelInput = {
    dx: number,
    dy: number,
    dz: number,
} & BaseInput

export type ClickInput = {
    x: number,
    y: number,
    pageX: number,
    pageY: number,
} & BaseInput

export type MoveInput = {
    x: number,
    y: number,
    pageX: number,
    pageY: number,
    inside: boolean,
} & BaseInput

export type PinchInput = {
    delta: number,
    distance: number,
    isStart: boolean
}

export type ResizeInput = {

}

const enum DraggingState {
    Stopped = 0,
    Started = 1,
    Moving = 2
}

type PointerEvent = {
    clientX: number
    clientY: number
    pageX: number
    pageY: number
}

interface InputObserver {
    noScroll: boolean
    noContextMenu: boolean

    drag: Subject<DragInput>,
    wheel: Subject<WheelInput>,
    pinch: Subject<PinchInput>,
    click: Subject<ClickInput>,
    move: Subject<MoveInput>,
    resize: Subject<ResizeInput>,

    dispose: () => void
}

namespace InputObserver {
    export function create (element: Element, props: InputObserverProps = {}): InputObserver {
        let { noScroll, noContextMenu } = { ...DefaultInputObserverProps, ...props }

        const lineHeight = toPixels('ex', element)

        let lastTouchDistance = 0
        const pointerStart = Vec2.zero()
        const pointerEnd = Vec2.zero()
        const pointerDelta = Vec2.zero()
        const rectSize = Vec2.zero()
        const modifiers: ModifiersKeys = {
            shift: false,
            alt: false,
            control: false,
            meta: false
        }

        let dragging: DraggingState = DraggingState.Stopped
        let disposed = false
        let buttons = 0

        const drag = new Subject<DragInput>()
        const click = new Subject<ClickInput>()
        const move = new Subject<MoveInput>()
        const wheel = new Subject<WheelInput>()
        const pinch = new Subject<PinchInput>()
        const resize = new Subject<ResizeInput>()

        attach()

        return {
            get noScroll () { return noScroll },
            set noScroll (value: boolean) { noScroll = value },
            get noContextMenu () { return noContextMenu },
            set noContextMenu (value: boolean) { noContextMenu = value },

            drag,
            wheel,
            pinch,
            click,
            move,
            resize,

            dispose
        }

        function attach () {
            element.addEventListener( 'contextmenu', onContextMenu, false )

            element.addEventListener('wheel', onMouseWheel, false)
            element.addEventListener('mousedown', onPointerDown as any, false)
            // for dragging to work outside canvas bounds,
            // mouse move/up events have to be added to a parent, i.e. window
            window.addEventListener('mousemove', onMouseMove as any, false)
            window.addEventListener('mouseup', onPointerUp as any, false)

            element.addEventListener('touchstart', onTouchStart as any, false)
            element.addEventListener('touchmove', onTouchMove as any, false)
            element.addEventListener('touchend', onTouchEnd as any, false)

            element.addEventListener('blur', handleBlur)
            element.addEventListener('keyup', handleMods as EventListener)
            element.addEventListener('keydown', handleMods as EventListener)
            element.addEventListener('keypress', handleMods as EventListener)

            window.addEventListener('resize', onResize, false)
        }

        function dispose () {
            if (disposed) return
            disposed = true

            element.removeEventListener( 'contextmenu', onContextMenu, false )

            element.removeEventListener('wheel', onMouseWheel, false)
            element.removeEventListener('mousedown', onMouseDown as any, false)
            window.removeEventListener('mousemove', onMouseMove as any, false)
            window.removeEventListener('mouseup', onMouseUp as any, false)

            element.removeEventListener('touchstart', onTouchStart as any, false)
            element.removeEventListener('touchmove', onTouchMove as any, false)
            element.removeEventListener('touchend', onTouchEnd as any, false)

            element.removeEventListener('blur', handleBlur)
            element.removeEventListener('keyup', handleMods as EventListener)
            element.removeEventListener('keydown', handleMods as EventListener)
            element.removeEventListener('keypress', handleMods as EventListener)

            window.removeEventListener('resize', onResize, false)
        }

        function preventDefault (ev: Event | Touch) {
            if ('preventDefault' in ev) ev.preventDefault()
        }

        function onContextMenu(event: Event) {
            if (noContextMenu) {
                event.preventDefault()
            }
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

        function getCenterTouch (ev: TouchEvent): PointerEvent {
            const t0 = ev.touches[0]
            const t1 = ev.touches[1]
            return {
                clientX: (t0.clientX + t1.clientX) / 2,
                clientY: (t0.clientY + t1.clientY) / 2,
                pageX: (t0.pageX + t1.pageX) / 2,
                pageY: (t0.pageY + t1.pageY) / 2
            }
        }

        function getTouchDistance (ev: TouchEvent) {
            const dx = ev.touches[0].pageX - ev.touches[1].pageX;
            const dy = ev.touches[0].pageY - ev.touches[1].pageY;
            return Math.sqrt(dx * dx + dy * dy);
        }

        function onTouchStart (ev: TouchEvent) {
            preventDefault(ev)

            if (ev.touches.length === 1) {
                buttons = ButtonsFlag.Primary
                onPointerDown(ev.touches[0])
            } else if (ev.touches.length >= 2) {
                buttons = ButtonsFlag.Secondary
                onPointerDown(getCenterTouch(ev))

                pinch.next({ distance: lastTouchDistance, delta: 0, isStart: true })
            }
        }

        function onTouchEnd (ev: TouchEvent) {
            preventDefault(ev)
        }

        function onTouchMove (ev: TouchEvent) {
            preventDefault(ev)

            if (ev.touches.length === 1) {
                buttons = ButtonsFlag.Primary
                onPointerMove(ev.touches[0])
            } else if (ev.touches.length >= 2) {
                const touchDistance = getTouchDistance(ev)
                if (lastTouchDistance - touchDistance < 4) {
                    buttons = ButtonsFlag.Secondary
                    onPointerMove(getCenterTouch(ev))
                } else {
                    pinch.next({
                        delta: lastTouchDistance - touchDistance,
                        distance: touchDistance,
                        isStart: false
                    })
                }
                lastTouchDistance = touchDistance
            }
        }

        function onMouseDown (ev: MouseEvent) {
            preventDefault(ev)

            buttons = getButtons(ev)
            onPointerDown(ev)
        }

        function onMouseMove (ev: MouseEvent) {
            preventDefault(ev)

            buttons = getButtons(ev)
            onPointerMove(ev)
        }

        function onMouseUp (ev: MouseEvent) {
            preventDefault(ev)

            buttons = getButtons(ev)
            onPointerUp(ev)
        }

        function onPointerDown (ev: PointerEvent) {
            eventOffset(pointerStart, ev)
            if (insideBounds(pointerStart)) {
                dragging = DraggingState.Started
            }
        }

        function onPointerUp (ev: PointerEvent) {
            dragging = DraggingState.Stopped

            if (Vec2.distance(pointerEnd, pointerStart) < 4) {
                eventOffset(pointerEnd, ev)

                const { pageX, pageY } = ev
                const [ x, y ] = pointerEnd

                click.next({ x, y, pageX, pageY, buttons, modifiers })
            }
        }

        function onPointerMove (ev: PointerEvent) {
            eventOffset(pointerEnd, ev)
            const { pageX, pageY } = ev
            const [ x, y ] = pointerEnd
            const inside = insideBounds(pointerEnd)
            move.next({ x, y, pageX, pageY, buttons, modifiers, inside })

            if (dragging === DraggingState.Stopped) return

            Vec2.div(pointerDelta, Vec2.sub(pointerDelta, pointerEnd, pointerStart), getClientSize(rectSize))

            const isStart = dragging === DraggingState.Started
            const [ dx, dy ] = pointerDelta
            drag.next({ x, y, dx, dy, pageX, pageY, buttons, modifiers, isStart })

            Vec2.copy(pointerStart, pointerEnd)
            dragging = DraggingState.Moving
        }

        function onMouseWheel(ev: MouseWheelEvent) {
            if (noScroll) {
                ev.preventDefault()
            }
            const mode = ev.deltaMode
            let dx = ev.deltaX || 0
            let dy = ev.deltaY || 0
            let dz = ev.deltaZ || 0
            let scale = 1
            switch (mode) {
                case 1: scale = lineHeight; break
                case 2: scale = window.innerHeight; break
            }
            scale *= 0.0001
            dx *= scale
            dy *= scale
            dz *= scale
            if (dx || dy || dz) {
                wheel.next({ dx, dy, dz, buttons, modifiers })
            }
        }

        function onResize (ev: Event) {
            resize.next()
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
            out[0] = element.clientWidth
            out[1] = element.clientHeight
            return out
        }

        function eventOffset (out: Vec2, ev: PointerEvent) {
            const cx = ev.clientX || 0
            const cy = ev.clientY || 0
            const rect = element.getBoundingClientRect()
            out[0] = cx - rect.left
            out[1] = cy - rect.top
            return out
        }
    }
}

export default InputObserver