/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified (use TypeScript, RxJS) from https://github.com/Jam3/touch-pinch,
 * copyright (c) 2014 Matt DesLauriers. MIT License
 */

import { Subject } from 'rxjs';

import { Vec2 } from 'mol-math/linear-algebra';
import { eventOffset } from './event-offset'

interface Finger {
    position: Vec2,
    touch?: Touch
}

function Finger (): Finger {
    return {
        position: Vec2.zero(),
        touch: undefined
    }
}

interface TouchPinch {
    pinching: boolean
    fingers: (Finger|undefined)[]
    indexOfTouch: (touch: Touch) => number

    start: Subject<number>
    end: Subject<void>
    place: Subject<{ newTouch?: Touch, oldTouch?: Touch}>
    change: Subject<{ currentDistance: number, lastDistance: number }>
    lift: Subject<{ removed: Touch, otherTouch?: Touch }>

    dispose: () => void
}

namespace TouchPinch {
    export function create (target: Element): TouchPinch {
        const fingers: (Finger|undefined)[] = []
        let activeCount = 0

        let lastDistance = 0
        let ended = false
        let disposed = false

        const start = new Subject<number>()
        const end = new Subject<void>()
        const place = new Subject<{ newTouch?: Touch, oldTouch?: Touch}>()
        const change = new Subject<{ currentDistance: number, lastDistance: number }>()
        const lift = new Subject<{ removed: Touch, otherTouch?: Touch }>()

        target.addEventListener('touchstart', onTouchStart as any, false)
        target.addEventListener('touchmove', onTouchMove as any, false)
        target.addEventListener('touchend', onTouchRemoved as any, false)
        target.addEventListener('touchcancel', onTouchRemoved as any, false)

        return {
            get pinching() { return activeCount === 2 },
            fingers,
            indexOfTouch,

            start,
            end,
            place,
            change,
            lift,

            dispose
        }

        function indexOfTouch (touch: Touch) {
            const id = touch.identifier
            for (let i = 0; i < fingers.length; i++) {
                const finger = fingers[i]
                if (finger && finger.touch && finger.touch.identifier === id) {
                    return i
                }
            }
            return -1
        }

        function dispose () {
            if (disposed) return
            disposed = true
            activeCount = 0
            fingers[0] = undefined
            fingers[1] = undefined
            lastDistance = 0
            ended = false
            target.removeEventListener('touchstart', onTouchStart as any, false)
            target.removeEventListener('touchmove', onTouchMove as any, false)
            target.removeEventListener('touchend', onTouchRemoved as any, false)
            target.removeEventListener('touchcancel', onTouchRemoved as any, false)
        }

        function onTouchStart (ev: TouchEvent) {
            for (let i = 0; i < ev.changedTouches.length; i++) {
                const newTouch = ev.changedTouches[i]
                const idx = indexOfTouch(newTouch)

                if (idx === -1 && activeCount < 2) {
                    const first = activeCount === 0

                    // newest and previous finger (previous may be undefined)
                    const newIndex = fingers[0] ? 1 : 0
                    const oldIndex = fingers[0] ? 0 : 1
                    const newFinger = Finger()

                    // add to stack
                    fingers[newIndex] = newFinger
                    activeCount++

                    // update touch event & position
                    newFinger.touch = newTouch
                    eventOffset(newFinger.position, newTouch, target)

                    const finger = fingers[oldIndex]
                    const oldTouch = finger ? finger.touch : undefined
                    place.next({ newTouch, oldTouch })

                    if (!first) {
                        const initialDistance = computeDistance()
                        ended = false
                        start.next(initialDistance)
                        lastDistance = initialDistance
                    }
                }
            }
        }

        function onTouchMove (ev: TouchEvent) {
            let changed = false
            for (let i = 0; i < ev.changedTouches.length; i++) {
                const movedTouch = ev.changedTouches[i]
                const idx = indexOfTouch(movedTouch)
                if (idx !== -1) {
                    const finger = fingers[idx]
                    if (finger) {
                        changed = true
                        finger.touch = movedTouch // avoid caching touches
                        eventOffset(finger.position, movedTouch, target)
                    }
                }
                }

                if (activeCount === 2 && changed) {
                const currentDistance = computeDistance()
                change.next({ currentDistance, lastDistance })
                lastDistance = currentDistance
            }
        }

        function onTouchRemoved (ev: TouchEvent) {
            for (let i = 0; i < ev.changedTouches.length; i++) {
                const removed = ev.changedTouches[i]
                const idx = indexOfTouch(removed)
                if (idx !== -1) {
                    fingers[idx] = undefined
                    activeCount--
                    const otherIdx = idx === 0 ? 1 : 0
                    const finger = fingers[otherIdx]
                    if (finger) {
                        const otherTouch = finger ? finger.touch : undefined
                        lift.next({ removed, otherTouch })
                    }
                }
            }

            if (!ended && activeCount !== 2) {
                ended = true
                end.next()
            }
        }

        function computeDistance () {
            const [ f1, f2 ] = fingers
            return (f1 && f2) ? Vec2.distance(f1.position, f2.position) : 0
        }
    }
}

export default TouchPinch