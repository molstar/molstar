/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Subject, Observable } from 'rxjs';

import { Vec2, EPSILON } from '../../mol-math/linear-algebra';

import { BitFlags, noop } from '../../mol-util';

export function getButtons(event: MouseEvent | Touch) {
    if (typeof event === 'object') {
        if ('buttons' in event) {
            return event.buttons;
        } else if ('which' in event) {
            const b = (event as any).which;  // 'any' to support older browsers
            if (b === 2) {
                return 4;
            } else if (b === 3) {
                return 2;
            } else if (b > 0) {
                return 1 << (b - 1);
            }
        }
    }
    return 0;
}

export function getButton(event: MouseEvent | Touch) {
    if (typeof event === 'object') {
        if ('button' in event) {
            const b = event.button;
            if (b === 1) {
                return 4;
            } else if (b === 2) {
                return 2;
            } else if (b >= 0) {
                return 1 << b;
            }
        }
    }
    return 0;
}

export function getModifiers(event: MouseEvent | Touch): ModifiersKeys {
    return {
        alt: 'altKey' in event ? event.altKey : false,
        shift: 'shiftKey' in event ? event.shiftKey : false,
        control: 'ctrlKey' in event ? event.ctrlKey : false,
        meta: 'metaKey' in event ? event.metaKey : false
    };
}

export const DefaultInputObserverProps = {
    noScroll: true,
    noMiddleClickScroll: true,
    noContextMenu: true,
    noPinchZoom: true
};
export type InputObserverProps = Partial<typeof DefaultInputObserverProps>

export type ModifiersKeys = {
    shift: boolean,
    alt: boolean,
    control: boolean,
    meta: boolean
}
export namespace ModifiersKeys {
    export const None = create();

    export function areEqual(a: ModifiersKeys, b: ModifiersKeys) {
        return a.shift === b.shift && a.alt === b.alt && a.control === b.control && a.meta === b.meta;
    }

    export function size(a?: ModifiersKeys) {
        if (!a) return 0;
        let ret = 0;
        if (!!a.shift) ret++;
        if (!!a.alt) ret++;
        if (!!a.control) ret++;
        if (!!a.meta) ret++;
        return ret;
    }

    export function create(modifierKeys: Partial<ModifiersKeys> = {}): ModifiersKeys {
        return {
            shift: !!modifierKeys.shift,
            alt: !!modifierKeys.alt,
            control: !!modifierKeys.control,
            meta: !!modifierKeys.meta
        };
    }
}

export type ButtonsType = BitFlags<ButtonsType.Flag>

export namespace ButtonsType {
    export const has: (btn: ButtonsType, f: Flag) => boolean = BitFlags.has;
    export const create: (fs: Flag) => ButtonsType = BitFlags.create;

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
    buttons: ButtonsType
    button: ButtonsType.Flag
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
    fraction: number,
    distance: number,
    isStart: boolean
} & BaseInput

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

    readonly drag: Observable<DragInput>,
    // Equivalent to mouseUp and touchEnd
    readonly interactionEnd: Observable<undefined>,
    readonly wheel: Observable<WheelInput>,
    readonly pinch: Observable<PinchInput>,
    readonly click: Observable<ClickInput>,
    readonly move: Observable<MoveInput>,
    readonly leave: Observable<undefined>,
    readonly enter: Observable<undefined>,
    readonly resize: Observable<ResizeInput>,
    readonly modifiers: Observable<ModifiersKeys>

    dispose: () => void
}

function createEvents() {
    return {
        drag: new Subject<DragInput>(),
        interactionEnd: new Subject<undefined>(),
        click: new Subject<ClickInput>(),
        move: new Subject<MoveInput>(),
        wheel: new Subject<WheelInput>(),
        pinch: new Subject<PinchInput>(),
        resize: new Subject<ResizeInput>(),
        leave: new Subject<undefined>(),
        enter: new Subject<undefined>(),
        modifiers: new Subject<ModifiersKeys>(),
    };
}

namespace InputObserver {
    export function create(props: InputObserverProps = {}): InputObserver {
        const { noScroll, noContextMenu } = { ...DefaultInputObserverProps, ...props };
        return {
            noScroll,
            noContextMenu,

            ...createEvents(),

            dispose: noop
        };
    }

    export function fromElement(element: Element, props: InputObserverProps = {}): InputObserver {
        let { noScroll, noMiddleClickScroll, noContextMenu, noPinchZoom } = { ...DefaultInputObserverProps, ...props };

        let lastTouchDistance = 0;
        const pointerDown = Vec2.zero();
        const pointerStart = Vec2.zero();
        const pointerEnd = Vec2.zero();
        const pointerDelta = Vec2.zero();
        const rectSize = Vec2.zero();
        const modifierKeys: ModifiersKeys = {
            shift: false,
            alt: false,
            control: false,
            meta: false
        };

        function getModifierKeys(): ModifiersKeys {
            return { ...modifierKeys };
        }

        let dragging: DraggingState = DraggingState.Stopped;
        let disposed = false;
        let buttons = ButtonsType.create(ButtonsType.Flag.None);
        let button = ButtonsType.Flag.None;
        let isInside = false;

        const events = createEvents();
        const { drag, interactionEnd, wheel, pinch, click, move, leave, enter, resize, modifiers } = events;

        attach();

        return {
            get noScroll () { return noScroll; },
            set noScroll (value: boolean) { noScroll = value; },
            get noContextMenu () { return noContextMenu; },
            set noContextMenu (value: boolean) { noContextMenu = value; },

            ...events,

            dispose
        };

        function attach() {
            element.addEventListener('contextmenu', onContextMenu, false );

            element.addEventListener('wheel', onMouseWheel as any, false);
            element.addEventListener('mousedown', onMouseDown as any, false);

            // for dragging to work outside canvas bounds,
            // mouse move/up events have to be added to a parent, i.e. window
            window.addEventListener('mousemove', onMouseMove as any, false);
            window.addEventListener('mouseup', onMouseUp as any, false);

            element.addEventListener('mouseenter', onMouseEnter as any, false);
            element.addEventListener('mouseleave', onMouseLeave as any, false);

            element.addEventListener('touchstart', onTouchStart as any, false);
            element.addEventListener('touchmove', onTouchMove as any, false);
            element.addEventListener('touchend', onTouchEnd as any, false);

            // reset buttons and modifier keys state when browser window looses focus
            window.addEventListener('blur', handleBlur);
            window.addEventListener('keyup', handleKeyUp as EventListener, false);
            window.addEventListener('keydown', handleKeyDown as EventListener, false);

            window.addEventListener('resize', onResize, false);
        }

        function dispose() {
            if (disposed) return;
            disposed = true;

            element.removeEventListener( 'contextmenu', onContextMenu, false );

            element.removeEventListener('wheel', onMouseWheel as any, false);
            element.removeEventListener('mousedown', onMouseDown as any, false);
            window.removeEventListener('mousemove', onMouseMove as any, false);
            window.removeEventListener('mouseup', onMouseUp as any, false);

            element.removeEventListener('mouseenter', onMouseEnter as any, false);
            element.removeEventListener('mouseleave', onMouseLeave as any, false);

            element.removeEventListener('touchstart', onTouchStart as any, false);
            element.removeEventListener('touchmove', onTouchMove as any, false);
            element.removeEventListener('touchend', onTouchEnd as any, false);

            window.removeEventListener('blur', handleBlur);
            window.removeEventListener('keyup', handleKeyUp as EventListener, false);
            window.removeEventListener('keydown', handleKeyDown as EventListener, false);

            window.removeEventListener('resize', onResize, false);
        }

        function onContextMenu(event: Event) {
            if (noContextMenu) {
                event.preventDefault();
            }
        }

        function updateModifierKeys(event: MouseEvent | WheelEvent | TouchEvent) {
            modifierKeys.alt = event.altKey;
            modifierKeys.shift = event.shiftKey;
            modifierKeys.control = event.ctrlKey;
            modifierKeys.meta = event.metaKey;
        }

        function handleBlur() {
            if (buttons || modifierKeys.shift || modifierKeys.alt || modifierKeys.meta || modifierKeys.control) {
                buttons = 0 as ButtonsType;
                modifierKeys.shift = modifierKeys.alt = modifierKeys.control = modifierKeys.meta = false;
            }
        }

        function handleKeyDown(event: KeyboardEvent) {
            let changed = false;
            if (!modifierKeys.alt && event.altKey) { changed = true; modifierKeys.alt = true; }
            if (!modifierKeys.shift && event.shiftKey) { changed = true; modifierKeys.shift = true; }
            if (!modifierKeys.control && event.ctrlKey) { changed = true; modifierKeys.control = true; }
            if (!modifierKeys.meta && event.metaKey) { changed = true; modifierKeys.meta = true; }

            if (changed && isInside) modifiers.next(getModifierKeys());
        }

        function handleKeyUp(event: KeyboardEvent) {
            let changed = false;

            if (modifierKeys.alt && !event.altKey) { changed = true; modifierKeys.alt = false; }
            if (modifierKeys.shift && !event.shiftKey) { changed = true; modifierKeys.shift = false; }
            if (modifierKeys.control && !event.ctrlKey) { changed = true; modifierKeys.control = false; }
            if (modifierKeys.meta && !event.metaKey) { changed = true; modifierKeys.meta = false; }

            if (changed && isInside) modifiers.next(getModifierKeys());
        }

        function getCenterTouch(ev: TouchEvent): PointerEvent {
            const t0 = ev.touches[0];
            const t1 = ev.touches[1];
            return {
                clientX: (t0.clientX + t1.clientX) / 2,
                clientY: (t0.clientY + t1.clientY) / 2,
                pageX: (t0.pageX + t1.pageX) / 2,
                pageY: (t0.pageY + t1.pageY) / 2
            };
        }

        function getTouchDistance(ev: TouchEvent) {
            const dx = ev.touches[0].pageX - ev.touches[1].pageX;
            const dy = ev.touches[0].pageY - ev.touches[1].pageY;
            return Math.sqrt(dx * dx + dy * dy);
        }

        function onTouchStart(ev: TouchEvent) {
            if (ev.touches.length === 1) {
                buttons = button = ButtonsType.Flag.Primary;
                onPointerDown(ev.touches[0]);
            } else if (ev.touches.length === 2) {
                buttons = ButtonsType.Flag.Secondary & ButtonsType.Flag.Auxilary;
                button = ButtonsType.Flag.Secondary;
                onPointerDown(getCenterTouch(ev));

                const touchDistance = getTouchDistance(ev);
                lastTouchDistance = touchDistance;
                pinch.next({
                    distance: touchDistance,
                    fraction: 1,
                    delta: 0,
                    isStart: true,
                    buttons,
                    button,
                    modifiers: getModifierKeys()
                });
            } else if (ev.touches.length === 3) {
                buttons = button = ButtonsType.Flag.Forth;
                onPointerDown(getCenterTouch(ev));
            }
        }

        function onTouchEnd(ev: TouchEvent) {
            endDrag();
        }

        function onTouchMove(ev: TouchEvent) {
            button = ButtonsType.Flag.None;

            if (noPinchZoom) {
                ev.preventDefault();
                ev.stopPropagation();
                if ((ev as any).originalEvent) {
                    (ev as any).originalEvent.preventDefault();
                    (ev as any).originalEvent.stopPropagation();
                }
            }

            if (ev.touches.length === 1) {
                buttons = ButtonsType.Flag.Primary;
                onPointerMove(ev.touches[0]);
            } else if (ev.touches.length === 2) {
                const touchDistance = getTouchDistance(ev);
                const touchDelta = lastTouchDistance - touchDistance;
                if (Math.abs(touchDelta) < 4) {
                    buttons = ButtonsType.Flag.Secondary;
                    onPointerMove(getCenterTouch(ev));
                } else {
                    buttons = ButtonsType.Flag.Auxilary;
                    updateModifierKeys(ev);
                    pinch.next({
                        delta: touchDelta,
                        fraction: lastTouchDistance / touchDistance,
                        distance: touchDistance,
                        isStart: false,
                        buttons,
                        button,
                        modifiers: getModifierKeys()
                    });
                }
                lastTouchDistance = touchDistance;
            } else if (ev.touches.length === 3) {
                buttons = ButtonsType.Flag.Forth;
                onPointerMove(getCenterTouch(ev));
            }
        }

        function onMouseDown(ev: MouseEvent) {
            updateModifierKeys(ev);
            buttons = getButtons(ev);
            button = getButton(ev);
            if (noMiddleClickScroll && buttons === ButtonsType.Flag.Auxilary) {
                ev.preventDefault;
            }
            onPointerDown(ev);
        }

        function onMouseMove(ev: MouseEvent) {
            updateModifierKeys(ev);
            buttons = getButtons(ev);
            button = ButtonsType.Flag.None;
            onPointerMove(ev);
        }

        function onMouseUp(ev: MouseEvent) {
            updateModifierKeys(ev);
            buttons = getButtons(ev);
            button = getButton(ev);
            onPointerUp(ev);
            endDrag();
        }

        function endDrag() {
            interactionEnd.next();
        }

        function onPointerDown(ev: PointerEvent) {
            eventOffset(pointerStart, ev);
            Vec2.copy(pointerDown, pointerStart);

            if (insideBounds(pointerStart)) {
                dragging = DraggingState.Started;
            }
        }

        function onPointerUp(ev: PointerEvent) {
            dragging = DraggingState.Stopped;

            eventOffset(pointerEnd, ev);
            if (Vec2.distance(pointerEnd, pointerDown) < 4) {
                const { pageX, pageY } = ev;
                const [ x, y ] = pointerEnd;

                click.next({ x, y, pageX, pageY, buttons, button, modifiers: getModifierKeys() });
            }
        }

        function onPointerMove(ev: PointerEvent) {
            eventOffset(pointerEnd, ev);
            const { pageX, pageY } = ev;
            const [ x, y ] = pointerEnd;
            const inside = insideBounds(pointerEnd);
            move.next({ x, y, pageX, pageY, buttons, button, modifiers: getModifierKeys(), inside });

            if (dragging === DraggingState.Stopped) return;

            Vec2.div(pointerDelta, Vec2.sub(pointerDelta, pointerEnd, pointerStart), getClientSize(rectSize));
            if (Vec2.magnitude(pointerDelta) < EPSILON) return;

            const isStart = dragging === DraggingState.Started;
            const [ dx, dy ] = pointerDelta;
            drag.next({ x, y, dx, dy, pageX, pageY, buttons, button, modifiers: getModifierKeys(), isStart });

            Vec2.copy(pointerStart, pointerEnd);
            dragging = DraggingState.Moving;
        }

        function onMouseWheel(ev: WheelEvent) {
            if (noScroll) {
                ev.preventDefault();
            }

            let scale = 1;
            switch (ev.deltaMode) {
                case 0: scale = 1; break; // pixels
                case 1: scale = 40; break; // lines
                case 2: scale = 800; break; // pages
            }

            const dx = (ev.deltaX || 0) * scale;
            const dy = (ev.deltaY || 0) * scale;
            const dz = (ev.deltaZ || 0) * scale;

            buttons = button = ButtonsType.Flag.Auxilary;

            if (dx || dy || dz) {
                wheel.next({ dx, dy, dz, buttons, button, modifiers: getModifierKeys() });
            }
        }

        function onMouseEnter(ev: Event) {
            isInside = true;
            enter.next();
        }

        function onMouseLeave(ev: Event) {
            isInside = false;
            leave.next();
        }

        function onResize(ev: Event) {
            resize.next();
        }

        function insideBounds(pos: Vec2) {
            if (element instanceof Window || element instanceof Document || element === document.body) {
                return true;
            } else {
                const rect = element.getBoundingClientRect();
                return pos[0] >= 0 && pos[1] >= 0 && pos[0] < rect.width && pos[1] < rect.height;
            }
        }

        function getClientSize(out: Vec2) {
            out[0] = element.clientWidth;
            out[1] = element.clientHeight;
            return out;
        }

        function eventOffset(out: Vec2, ev: PointerEvent) {
            const cx = ev.clientX || 0;
            const cy = ev.clientY || 0;
            const rect = element.getBoundingClientRect();
            out[0] = cx - rect.left;
            out[1] = cy - rect.top;
            return out;
        }
    }
}

export default InputObserver;