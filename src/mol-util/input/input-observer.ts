/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Subject } from 'rxjs';

import { Vec2 } from 'mol-math/linear-algebra';

import { BitFlags } from 'mol-util';

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
                return 1 << (b - 1)
            }
        } else if ('button' in event) {
            const b = (event as any).button  // 'any' to support older browsers
            if (b === 1) {
                return 4
            } else if (b === 2) {
                return 2
            } else if (b >= 0) {
                return 1 << b
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

export interface ButtonsType extends BitFlags<ButtonsType.Flag> { }

export namespace ButtonsType {
    export const has: (ss: ButtonsType, f: Flag) => boolean = BitFlags.has
    export const create: (fs: Flag) => ButtonsType = BitFlags.create

    export const enum Flag {
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
    // Equivalent to mouseUp and touchEnd
    interactionEnd: Subject<undefined>,
    wheel: Subject<WheelInput>,
    pinch: Subject<PinchInput>,
    click: Subject<ClickInput>,
    move: Subject<MoveInput>,
    leave: Subject<undefined>,
    enter: Subject<undefined>,
    resize: Subject<ResizeInput>,

    dispose: () => void
}

namespace InputObserver {
    export function create (element: Element, props: InputObserverProps = {}): InputObserver {
        let { noScroll, noContextMenu } = { ...DefaultInputObserverProps, ...props }

        let lastTouchDistance = 0
        const pointerDown = Vec2.zero()
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
        const interactionEnd = new Subject<undefined>();
        const click = new Subject<ClickInput>()
        const move = new Subject<MoveInput>()
        const wheel = new Subject<WheelInput>()
        const pinch = new Subject<PinchInput>()
        const resize = new Subject<ResizeInput>()
        const leave = new Subject<undefined>()
        const enter = new Subject<undefined>()

        attach()

        return {
            get noScroll () { return noScroll },
            set noScroll (value: boolean) { noScroll = value },
            get noContextMenu () { return noContextMenu },
            set noContextMenu (value: boolean) { noContextMenu = value },

            drag,
            interactionEnd,
            wheel,
            pinch,
            click,
            move,
            leave,
            enter,
            resize,

            dispose
        }

        function attach () {
            element.addEventListener( 'contextmenu', onContextMenu, false )

            element.addEventListener('wheel', onMouseWheel as any, false)
            element.addEventListener('mousedown', onMouseDown as any, false)
            // for dragging to work outside canvas bounds,
            // mouse move/up events have to be added to a parent, i.e. window
            window.addEventListener('mousemove', onMouseMove as any, false)
            window.addEventListener('mouseup', onMouseUp as any, false)

            element.addEventListener('mouseenter', onMouseEnter as any, false)
            element.addEventListener('mouseleave', onMouseLeave as any, false)

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

            element.removeEventListener('wheel', onMouseWheel as any, false)
            element.removeEventListener('mousedown', onMouseDown as any, false)
            window.removeEventListener('mousemove', onMouseMove as any, false)
            window.removeEventListener('mouseup', onMouseUp as any, false)

            element.removeEventListener('mouseenter', onMouseEnter as any, false)
            element.removeEventListener('mouseleave', onMouseLeave as any, false)

            element.removeEventListener('touchstart', onTouchStart as any, false)
            element.removeEventListener('touchmove', onTouchMove as any, false)
            element.removeEventListener('touchend', onTouchEnd as any, false)

            element.removeEventListener('blur', handleBlur)
            element.removeEventListener('keyup', handleMods as EventListener)
            element.removeEventListener('keydown', handleMods as EventListener)
            element.removeEventListener('keypress', handleMods as EventListener)

            window.removeEventListener('resize', onResize, false)
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
            if (ev.touches.length === 1) {
                buttons = ButtonsType.Flag.Primary
                onPointerDown(ev.touches[0])
            } else if (ev.touches.length >= 2) {
                buttons = ButtonsType.Flag.Secondary
                onPointerDown(getCenterTouch(ev))

                pinch.next({ distance: lastTouchDistance, delta: 0, isStart: true })
            }
        }

        function onTouchEnd (ev: TouchEvent) {
            endDrag()
        }

        function onTouchMove (ev: TouchEvent) {
            if (ev.touches.length === 1) {
                buttons = ButtonsType.Flag.Primary
                onPointerMove(ev.touches[0])
            } else if (ev.touches.length >= 2) {
                const touchDistance = getTouchDistance(ev)
                if (lastTouchDistance - touchDistance < 4) {
                    buttons = ButtonsType.Flag.Secondary
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
            buttons = getButtons(ev)
            onPointerDown(ev)
        }

        function onMouseMove (ev: MouseEvent) {
            buttons = getButtons(ev)
            onPointerMove(ev)
        }

        function onMouseUp (ev: MouseEvent) {
            onPointerUp(ev)
            endDrag()
        }

        function endDrag() {
            interactionEnd.next()
        }

        function onPointerDown (ev: PointerEvent) {
            eventOffset(pointerStart, ev)
            Vec2.copy(pointerDown, pointerStart)

            if (insideBounds(pointerStart)) {
                dragging = DraggingState.Started
            }
        }

        function onPointerUp (ev: PointerEvent) {
            dragging = DraggingState.Stopped

            eventOffset(pointerEnd, ev);
            if (Vec2.distance(pointerEnd, pointerDown) < 4) {
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

        function onMouseWheel(ev: WheelEvent) {
            if (noScroll) {
                ev.preventDefault()
            }

            let scale = 1
            switch (ev.deltaMode) {
                case 0: scale = 1; break // pixels
                case 1: scale = 40; break // lines
                case 2: scale = 800; break // pages
            }

            const dx = (ev.deltaX || 0) * scale
            const dy = (ev.deltaY || 0) * scale
            const dz = (ev.deltaZ || 0) * scale

            if (dx || dy || dz) {
                wheel.next({ dx, dy, dz, buttons, modifiers })
            }
        }

        function onMouseEnter (ev: Event) {
            enter.next();
        }

        function onMouseLeave (ev: Event) {
            leave.next();
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