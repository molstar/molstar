/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Russell Parker <russell@benchling.com>
 * @author Herman Bergwerf <post@hbergwerf.nl>
 */

import { Subject, Observable } from 'rxjs';
import { Viewport } from '../../mol-canvas3d/camera/util';
import { Vec2, EPSILON } from '../../mol-math/linear-algebra';
import { BitFlags, noop } from '../../mol-util';
import { Ray3D } from '../../mol-math/geometry/primitives/ray3d';

export function getButtons(event: MouseEvent | Touch) {
    if (typeof event === 'object') {
        if ('buttons' in event) {
            return event.buttons;
        } else if ('which' in event) {
            const b = (event as any).which; // 'any' to support older browsers
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
    noPinchZoom: true,
    noTextSelect: true,
    preventGestures: false,
    mask: (x: number, y: number) => true,

    pixelScale: 1
};
export type InputObserverProps = Partial<typeof DefaultInputObserverProps>

export type ModifiersKeys = {
    shift: boolean,
    alt: boolean,
    control: boolean,
    meta: boolean
}
export namespace ModifiersKeys {
    export const None: Readonly<ModifiersKeys> = create();

    export function areEqual(a: ModifiersKeys, b: ModifiersKeys) {
        return a.shift === b.shift && a.alt === b.alt && a.control === b.control && a.meta === b.meta;
    }

    export function areNone(a: ModifiersKeys) {
        return areEqual(a, None);
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

    export enum Flag {
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

export type KeyCode = string

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
    useDelta?: boolean
} & BaseInput

export type WheelInput = {
    x: number,
    y: number,
    pageX: number,
    pageY: number,
    dx: number,
    dy: number,
    dz: number,
    spinX: number,
    spinY: number
} & BaseInput

export type ClickInput = {
    x: number,
    y: number,
    pageX: number,
    pageY: number,
    ray?: Ray3D,
} & BaseInput

export type MoveInput = {
    x: number,
    y: number,
    pageX: number,
    pageY: number,
    ray?: Ray3D,
    movementX?: number,
    movementY?: number,
    inside: boolean,
    // Move is subscribed to window element
    // This indicates that the event originated from the element the InputObserver was created on
    onElement: boolean
} & BaseInput

export type PinchInput = {
    isStart: boolean,
    distance: number,
    delta: number,
    fraction: number,
    fractionDelta: number,
    startX: number,
    startY: number,
    centerPageX: number,
    centerPageY: number,
} & BaseInput

export type GestureInput = {
    scale: number,
    rotation: number,
    deltaScale: number,
    deltaRotation: number
    isStart?: boolean,
    isEnd?: boolean
}

export type KeyInput = {
    key: string,
    code: string,
    modifiers: ModifiersKeys
    x: number,
    y: number,
    pageX: number,
    pageY: number,
    /** for overwriting browser shortcuts like `ctrl+s` as needed */
    preventDefault: () => void
}

export const EmptyKeyInput: KeyInput = {
    key: '',
    code: '',
    modifiers: ModifiersKeys.None,
    x: -1,
    y: -1,
    pageX: -1,
    pageY: -1,
    preventDefault: noop,
};

export type ResizeInput = {

}

enum DraggingState {
    Stopped = 0,
    Started = 1,
    Moving = 2
}

type PointerEvent = {
    clientX: number
    clientY: number
    pageX: number
    pageY: number
    movementX?: number
    movementY?: number
    target: EventTarget | null

    preventDefault?: () => void
}

type GestureEvent = {
    scale: number,
    rotation: number,
} & MouseEvent

interface InputObserver {
    noScroll: boolean
    noContextMenu: boolean

    readonly width: number
    readonly height: number
    readonly pixelRatio: number
    readonly pointerLock: boolean

    readonly drag: Observable<DragInput>,
    // Equivalent to mouseUp and touchEnd
    readonly interactionEnd: Observable<undefined>,
    readonly wheel: Observable<WheelInput>,
    readonly pinch: Observable<PinchInput>,
    readonly gesture: Observable<GestureInput>,
    readonly click: Observable<ClickInput>,
    readonly move: Observable<MoveInput>,
    readonly leave: Observable<undefined>,
    readonly enter: Observable<undefined>,
    readonly resize: Observable<ResizeInput>,
    readonly modifiers: Observable<ModifiersKeys>
    readonly key: Observable<KeyInput>
    readonly keyUp: Observable<KeyInput>
    readonly keyDown: Observable<KeyInput>
    readonly lock: Observable<boolean>

    setPixelScale: (pixelScale: number) => void

    requestPointerLock: (viewport: Viewport) => void
    exitPointerLock: () => void
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
        gesture: new Subject<GestureInput>(),
        resize: new Subject<ResizeInput>(),
        leave: new Subject<undefined>(),
        enter: new Subject<undefined>(),
        modifiers: new Subject<ModifiersKeys>(),
        key: new Subject<KeyInput>(),
        keyUp: new Subject<KeyInput>(),
        keyDown: new Subject<KeyInput>(),
        lock: new Subject<boolean>(),
    };
}

const AllowedNonPrintableKeys = ['Backspace', 'Delete'];

namespace InputObserver {
    export function create(props: InputObserverProps = {}): InputObserver {
        const { noScroll, noContextMenu } = { ...DefaultInputObserverProps, ...props };
        return {
            noScroll,
            noContextMenu,
            pointerLock: false,

            width: 0,
            height: 0,
            pixelRatio: 1,

            ...createEvents(),

            setPixelScale: noop,

            requestPointerLock: noop,
            exitPointerLock: noop,
            dispose: noop
        };
    }

    export function fromElement(element: Element, props: InputObserverProps = {}): InputObserver {
        let { noScroll, noMiddleClickScroll, noContextMenu, noPinchZoom, noTextSelect, mask, pixelScale, preventGestures } = { ...DefaultInputObserverProps, ...props };

        let width = element.clientWidth * pixelRatio();
        let height = element.clientHeight * pixelRatio();

        let isLocked = false;
        let lockedViewport = Viewport();

        const pointerDown = Vec2();
        const pointerStart = Vec2();
        const pointerEnd = Vec2();
        const pointerDelta = Vec2();
        const rectSize = Vec2();
        const modifierKeys: ModifiersKeys = {
            shift: false,
            alt: false,
            control: false,
            meta: false
        };
        const position = {
            x: -1,
            y: -1,
            pageX: -1,
            pageY: -1,
        };

        function pixelRatio() {
            return window.devicePixelRatio * pixelScale;
        }

        function getModifierKeys(): ModifiersKeys {
            return { ...modifierKeys };
        }

        function getKeyOnElement(event: Event): boolean {
            return event.target === document.body || event.target === element;
        }

        let dragging: DraggingState = DraggingState.Stopped;
        let disposed = false;
        let buttons = ButtonsType.create(ButtonsType.Flag.None);
        let button = ButtonsType.Flag.None;
        let isInside = false;
        let hasMoved = false;

        let resizeObserver: ResizeObserver | undefined;
        if (typeof window.ResizeObserver !== 'undefined') {
            resizeObserver = new window.ResizeObserver(onResize);
        }

        const events = createEvents();
        const { drag, interactionEnd, wheel, pinch, gesture, click, move, leave, enter, resize, modifiers, key, keyUp, keyDown, lock } = events;

        attach();

        function attach() {
            element.addEventListener('contextmenu', onContextMenu as any, false);

            element.addEventListener('wheel', onMouseWheel as any, false);
            element.addEventListener('mousedown', onMouseDown as any, false);

            // for dragging to work outside canvas bounds,
            // mouse move/up events have to be added to a parent, i.e. window
            window.addEventListener('mousemove', onMouseMove as any, false);
            window.addEventListener('mouseup', onMouseUp as any, false);

            element.addEventListener('touchstart', onTouchStart as any, false);
            element.addEventListener('touchmove', onTouchMove as any, false);
            element.addEventListener('touchend', onTouchEnd as any, false);

            element.addEventListener('gesturechange', onGestureChange as any, false);
            element.addEventListener('gesturestart', onGestureStart as any, false);
            element.addEventListener('gestureend', onGestureEnd as any, false);

            // reset buttons and modifier keys state when browser window looses focus
            window.addEventListener('blur', handleBlur);
            window.addEventListener('keyup', handleKeyUp as EventListener, false);
            window.addEventListener('keydown', handleKeyDown as EventListener, false);
            window.addEventListener('keypress', handleKeyPress as EventListener, false);

            document.addEventListener('pointerlockchange', onPointerLockChange, false);
            document.addEventListener('pointerlockerror', onPointerLockError, false);

            if (resizeObserver != null) {
                resizeObserver.observe(element.parentElement!);
            } else {
                window.addEventListener('resize', onResize, false);
            }
        }

        function dispose() {
            if (disposed) return;
            disposed = true;

            element.removeEventListener('contextmenu', onContextMenu as any, false);

            element.removeEventListener('wheel', onMouseWheel as any, false);
            element.removeEventListener('mousedown', onMouseDown as any, false);
            window.removeEventListener('mousemove', onMouseMove as any, false);
            window.removeEventListener('mouseup', onMouseUp as any, false);

            element.removeEventListener('touchstart', onTouchStart as any, false);
            element.removeEventListener('touchmove', onTouchMove as any, false);
            element.removeEventListener('touchend', onTouchEnd as any, false);

            element.removeEventListener('gesturechange', onGestureChange as any, false);
            element.removeEventListener('gesturestart', onGestureStart as any, false);
            element.removeEventListener('gestureend', onGestureEnd as any, false);

            window.removeEventListener('blur', handleBlur);
            window.removeEventListener('keyup', handleKeyUp as EventListener, false);
            window.removeEventListener('keydown', handleKeyDown as EventListener, false);
            window.removeEventListener('keypress', handleKeyPress as EventListener, false);

            document.removeEventListener('pointerlockchange', onPointerLockChange, false);
            document.removeEventListener('pointerlockerror', onPointerLockError, false);

            cross.remove();

            if (resizeObserver != null) {
                resizeObserver.unobserve(element.parentElement!);
                resizeObserver.disconnect();
            } else {
                window.removeEventListener('resize', onResize, false);
            }
        }

        function onPointerLockChange() {
            if (element.ownerDocument.pointerLockElement === element) {
                isLocked = true;
            } else {
                isLocked = false;
            }
            toggleCross(isLocked);
            lock.next(isLocked);
        }

        function onPointerLockError() {
            console.error('Unable to use Pointer Lock API');
            isLocked = false;
            toggleCross(isLocked);
            lock.next(isLocked);
        }

        function onContextMenu(event: MouseEvent) {
            if (!mask(event.clientX, event.clientY)) return;

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

            if (getKeyOnElement(event) && isInside) {
                keyDown.next({
                    key: event.key,
                    code: event.code,
                    modifiers: getModifierKeys(),
                    ...position,
                    preventDefault: () => event.preventDefault(),
                });
            }
        }

        function handleKeyUp(event: KeyboardEvent) {
            let changed = false;

            if (modifierKeys.alt && !event.altKey) { changed = true; modifierKeys.alt = false; }
            if (modifierKeys.shift && !event.shiftKey) { changed = true; modifierKeys.shift = false; }
            if (modifierKeys.control && !event.ctrlKey) { changed = true; modifierKeys.control = false; }
            if (modifierKeys.meta && !event.metaKey) { changed = true; modifierKeys.meta = false; }

            if (changed && isInside) modifiers.next(getModifierKeys());

            if (AllowedNonPrintableKeys.includes(event.key)) handleKeyPress(event);

            if (getKeyOnElement(event) && isInside) {
                keyUp.next({
                    key: event.key,
                    code: event.code,
                    modifiers: getModifierKeys(),
                    ...position,
                    preventDefault: () => event.preventDefault(),
                });
            }
        }

        function handleKeyPress(event: KeyboardEvent) {
            if (!getKeyOnElement(event) || !isInside) return;

            key.next({
                key: event.key,
                code: event.code,
                modifiers: getModifierKeys(),
                ...position,
                preventDefault: () => event.preventDefault(),
            });
        }

        function getCenterTouch(ev: TouchEvent): PointerEvent {
            const t0 = ev.touches[0];
            const t1 = ev.touches[1];
            return {
                clientX: (t0.clientX + t1.clientX) / 2,
                clientY: (t0.clientY + t1.clientY) / 2,
                pageX: (t0.pageX + t1.pageX) / 2,
                pageY: (t0.pageY + t1.pageY) / 2,
                target: ev.target
            };
        }

        function getTouchDistance(ev: TouchEvent) {
            const dx = ev.touches[0].pageX - ev.touches[1].pageX;
            const dy = ev.touches[0].pageY - ev.touches[1].pageY;
            return Math.sqrt(dx * dx + dy * dy);
        }

        let singleTouchDistance = -1;
        let lastSingleTouch: Touch | undefined = undefined;
        const singleTouchPosition = Vec2(), singleTouchTmp = Vec2();

        function updateSingleTouchDistance(ev: TouchEvent) {
            if (singleTouchDistance < 0) return;

            Vec2.set(singleTouchTmp, ev.touches[0].pageX, ev.touches[0].pageY);
            singleTouchDistance += Vec2.distance(singleTouchPosition, singleTouchTmp);
            Vec2.copy(singleTouchPosition, singleTouchTmp);
        }

        const firstTouchStart = Vec2();
        let firstTouchStartSet = false;
        let initialTouchDistance = 0, lastTouchFraction = 1;

        function onTouchStart(ev: TouchEvent) {
            ev.preventDefault();

            lastSingleTouch = undefined;
            singleTouchDistance = -1;
            if (ev.touches.length === 1) {
                buttons = button = ButtonsType.Flag.Primary;

                singleTouchDistance = 0;
                Vec2.set(singleTouchPosition, ev.touches[0].pageX, ev.touches[0].pageY);
                lastSingleTouch = ev.touches[0];

                onPointerDown(ev.touches[0]);
                Vec2.copy(firstTouchStart, pointerStart);
                firstTouchStartSet = true;
            } else if (ev.touches.length === 2) {
                buttons = ButtonsType.Flag.Secondary | ButtonsType.Flag.Auxilary;
                button = ButtonsType.Flag.Secondary;
                updateModifierKeys(ev);

                lastTouchFraction = 1;
                initialTouchDistance = getTouchDistance(ev);
                const { pageX: centerPageX, pageY: centerPageY } = getPagePosition(getCenterTouch(ev));
                if (!firstTouchStartSet) {
                    eventOffset(firstTouchStart, getCenterTouch(ev));
                    firstTouchStartSet = true;
                }

                pinch.next({
                    isStart: true,
                    distance: initialTouchDistance,
                    delta: 0,
                    fraction: lastTouchFraction,
                    fractionDelta: 0,
                    startX: firstTouchStart[0],
                    startY: firstTouchStart[1],
                    centerPageX,
                    centerPageY,
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

            if (lastSingleTouch && singleTouchDistance <= 4) {
                const t = lastSingleTouch;
                if (!mask(t.clientX, t.clientY)) return;

                eventOffset(singleTouchTmp, t);
                const { pageX, pageY } = getPagePosition(t);
                const [x, y] = singleTouchTmp;

                click.next({ x, y, pageX, pageY, buttons, button, modifiers: getModifierKeys() });
            }
            lastSingleTouch = undefined;
            firstTouchStartSet = false;
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

            lastSingleTouch = undefined;
            if (ev.touches.length === 1) {
                buttons = ButtonsType.Flag.Primary;
                lastSingleTouch = ev.touches[0];
                updateSingleTouchDistance(ev);
                onPointerMove(ev.touches[0]);
            } else if (ev.touches.length === 2) {
                buttons = ButtonsType.Flag.Secondary | ButtonsType.Flag.Auxilary;
                button = ButtonsType.Flag.Secondary;
                updateModifierKeys(ev);

                const { pageX: centerPageX, pageY: centerPageY } = getPagePosition(getCenterTouch(ev));
                const distance = getTouchDistance(ev);
                const delta = initialTouchDistance - distance;
                const fraction = initialTouchDistance / distance;
                const fractionDelta = fraction - lastTouchFraction;
                lastTouchFraction = fraction;

                pinch.next({
                    isStart: false,
                    distance,
                    delta,
                    fraction,
                    fractionDelta,
                    startX: firstTouchStart[0],
                    startY: firstTouchStart[1],
                    centerPageX,
                    centerPageY,
                    buttons,
                    button,
                    modifiers: getModifierKeys()
                });
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
            interactionEnd.next(void 0);
        }

        function onPointerDown(ev: PointerEvent) {
            if (!mask(ev.clientX, ev.clientY)) return;

            eventOffset(pointerStart, ev);
            Vec2.copy(pointerDown, pointerStart);

            if (insideBounds(pointerStart)) {
                dragging = DraggingState.Started;
            }
        }

        function onPointerUp(ev: PointerEvent) {
            dragging = DraggingState.Stopped;
            if (!mask(ev.clientX, ev.clientY)) return;

            eventOffset(pointerEnd, ev);
            if (!hasMoved && Vec2.distance(pointerEnd, pointerDown) < 4) {
                const { pageX, pageY } = getPagePosition(ev);
                const [x, y] = pointerEnd;

                click.next({ x, y, pageX, pageY, buttons, button, modifiers: getModifierKeys() });
            }
            hasMoved = false;
        }

        function onPointerMove(ev: PointerEvent) {
            eventOffset(pointerEnd, ev);
            const { pageX, pageY } = getPagePosition(ev);
            const [x, y] = pointerEnd;
            const { movementX, movementY } = ev;

            const inside = insideBounds(pointerEnd) && mask(ev.clientX, ev.clientY);
            if (isInside && !inside) {
                leave.next(void 0);
            } else if (!isInside && inside) {
                enter.next(void 0);
            }
            isInside = inside;

            position.x = x;
            position.y = y;
            position.pageX = pageX;
            position.pageY = pageY;

            move.next({ x, y, pageX, pageY, movementX, movementY, buttons, button, modifiers: getModifierKeys(), inside, onElement: ev.target === element });

            if (dragging === DraggingState.Stopped) return;

            if (noTextSelect) {
                ev.preventDefault?.();
            }

            Vec2.div(pointerDelta, Vec2.sub(pointerDelta, pointerEnd, pointerStart), getClientSize(rectSize));
            if (Vec2.magnitude(pointerDelta) < EPSILON) return;

            const isStart = dragging === DraggingState.Started;
            if (isStart && !mask(ev.clientX, ev.clientY)) return;

            if (Vec2.distance(pointerEnd, pointerDown) >= 4) {
                hasMoved = true;
            }

            const [dx, dy] = pointerDelta;
            drag.next({ x, y, dx, dy, pageX, pageY, buttons, button, modifiers: getModifierKeys(), isStart });

            Vec2.copy(pointerStart, pointerEnd);
            dragging = DraggingState.Moving;
        }

        function onMouseWheel(ev: WheelEvent) {
            if (!mask(ev.clientX, ev.clientY)) return;

            eventOffset(pointerEnd, ev);
            const { pageX, pageY } = getPagePosition(ev);
            const [x, y] = pointerEnd;

            if (noScroll) {
                ev.preventDefault();
            }

            const normalized = normalizeWheel(ev);
            buttons = button = ButtonsType.Flag.Auxilary;

            if (normalized.dx || normalized.dy || normalized.dz) {
                wheel.next({ x, y, pageX, pageY, ...normalized, buttons, button, modifiers: getModifierKeys() });
            }
        }

        function tryPreventGesture(ev: GestureEvent) {
            // console.log(ev, preventGestures);
            if (!preventGestures) return;
            ev.preventDefault();
            ev.stopImmediatePropagation?.();
            ev.stopPropagation?.();
        }

        let prevGestureScale = 0, prevGestureRotation = 0;

        function onGestureStart(ev: GestureEvent) {
            tryPreventGesture(ev);
            prevGestureScale = ev.scale;
            prevGestureRotation = ev.rotation;
            gesture.next({ scale: ev.scale, rotation: ev.rotation, deltaRotation: 0, deltaScale: 0, isStart: true });
        }

        function gestureDelta(ev: GestureEvent, isEnd?: boolean) {
            gesture.next({
                scale: ev.scale,
                rotation: ev.rotation,
                deltaRotation: prevGestureRotation - ev.rotation,
                deltaScale: prevGestureScale - ev.scale,
                isEnd
            });
            prevGestureRotation = ev.rotation;
            prevGestureScale = ev.scale;
        }

        function onGestureChange(ev: GestureEvent) {
            tryPreventGesture(ev);
            gestureDelta(ev);
        }

        function onGestureEnd(ev: GestureEvent) {
            tryPreventGesture(ev);
            gestureDelta(ev, true);
        }

        function onResize() {
            width = element.clientWidth * pixelRatio();
            height = element.clientHeight * pixelRatio();
            resize.next({});
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

        function eventOffset(out: Vec2, ev: { clientX: number, clientY: number }) {
            width = element.clientWidth * pixelRatio();
            height = element.clientHeight * pixelRatio();

            if (isLocked) {
                const pr = pixelRatio();
                out[0] = (lockedViewport.x + lockedViewport.width / 2) / pr;
                out[1] = (height - (lockedViewport.y + lockedViewport.height / 2)) / pr;
            } else {
                const rect = element.getBoundingClientRect();
                out[0] = (ev.clientX || 0) - rect.left;
                out[1] = (ev.clientY || 0) - rect.top;
            }
            return out;
        }

        function getPagePosition(ev: { pageX: number, pageY: number }) {
            if (isLocked) {
                return {
                    pageX: Math.round(window.innerWidth / 2) + lockedViewport.x,
                    pageY: Math.round(window.innerHeight / 2) + lockedViewport.y
                };
            } else {
                return {
                    pageX: ev.pageX,
                    pageY: ev.pageY
                };
            }
        }

        const crossWidth = 30;
        const cross = addCross();

        function addCross() {
            const cross = document.createElement('div');

            const b = '30%';
            const t = '10%';
            const c = `#000 ${b}, #0000 0 calc(100% - ${b}), #000 0`;
            const vline = `linear-gradient(0deg, ${c}) 50%/${t} 100% no-repeat`;
            const hline = `linear-gradient(90deg, ${c}) 50%/100% ${t} no-repeat`;
            const cdot = 'radial-gradient(circle at 50%, #000 5%, #0000 5%)';
            Object.assign(cross.style, {
                width: `${crossWidth}px`,
                aspectRatio: 1,
                background: `${vline}, ${hline}, ${cdot}`,
                display: 'none',
                zIndex: 1000,
                position: 'absolute',
                mixBlendMode: 'difference',
                filter: 'invert(1)',
            });

            element.parentElement?.appendChild(cross);

            return cross;
        }

        function toggleCross(value: boolean) {
            cross.style.display = value ? 'block' : 'none';
            if (value) {
                const pr = pixelRatio();
                const offsetX = (lockedViewport.x + lockedViewport.width / 2) / pr;
                const offsetY = (lockedViewport.y + lockedViewport.height / 2) / pr;
                cross.style.width = `${crossWidth}px`;
                cross.style.left = `calc(${offsetX}px - ${crossWidth / 2}px)`;
                cross.style.bottom = `calc(${offsetY}px - ${crossWidth / 2}px)`;
            }
        }

        return {
            get noScroll() { return noScroll; },
            set noScroll(value: boolean) { noScroll = value; },
            get noContextMenu() { return noContextMenu; },
            set noContextMenu(value: boolean) { noContextMenu = value; },

            get width() { return width; },
            get height() { return height; },
            get pixelRatio() { return pixelRatio(); },
            get pointerLock() { return isLocked; },

            ...events,

            setPixelScale: (value: number) => {
                pixelScale = value;
                width = element.clientWidth * pixelRatio();
                height = element.clientHeight * pixelRatio();
            },

            requestPointerLock: (viewport: Viewport) => {
                lockedViewport = viewport;
                if (!isLocked) {
                    element.requestPointerLock();
                }
            },
            exitPointerLock: () => {
                if (isLocked) {
                    element.ownerDocument.exitPointerLock();
                }
            },
            dispose
        };
    }
}


// Adapted from https://stackoverflow.com/a/30134826
// License: https://creativecommons.org/licenses/by-sa/3.0/
export function normalizeWheel(event: any) {
    // Reasonable defaults
    const PIXEL_STEP = 10;
    const LINE_HEIGHT = 40;
    const PAGE_HEIGHT = 800;
    let spinX = 0, spinY = 0,
        dx = 0, dy = 0, dz = 0; // pixelX, pixelY, pixelZ

    // Legacy
    if ('detail' in event) { spinY = event.detail; }
    if ('wheelDelta' in event) { spinY = -event.wheelDelta / 120; }
    if ('wheelDeltaY' in event) { spinY = -event.wheelDeltaY / 120; }
    if ('wheelDeltaX' in event) { spinX = -event.wheelDeltaX / 120; }

    // side scrolling on FF with DOMMouseScroll
    if ('axis' in event && event.axis === event.HORIZONTAL_AXIS) {
        spinX = spinY;
        spinY = 0;
    }

    dx = spinX * PIXEL_STEP;
    dy = spinY * PIXEL_STEP;

    if ('deltaY' in event) { dy = event.deltaY; }
    if ('deltaX' in event) { dx = event.deltaX; }
    if ('deltaZ' in event) { dz = event.deltaZ; }

    if ((dx || dy || dz) && event.deltaMode) {
        if (event.deltaMode === 1) { // delta in LINE units
            dx *= LINE_HEIGHT;
            dy *= LINE_HEIGHT;
            dz *= LINE_HEIGHT;
        } else { // delta in PAGE units
            dx *= PAGE_HEIGHT;
            dy *= PAGE_HEIGHT;
            dz *= PAGE_HEIGHT;
        }
    }

    // Fall-back if spin cannot be determined
    if (dx && !spinX) { spinX = (dx < 1) ? -1 : 1; }
    if (dy && !spinY) { spinY = (dy < 1) ? -1 : 1; }

    return { spinX, spinY, dx, dy, dz };
}


export { InputObserver };