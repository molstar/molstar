/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/mikolalysenko/mouse-event,
 * copyright (c) 2015 Mikola Lysenko. MIT License
 */

export function buttons(event: MouseEvent) {
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

export function element(event: MouseEvent) {
    return event.target as Element
}

export function x(event: MouseEvent) {
    if (typeof event === 'object') {
        if ('offsetX' in event) {
            return event.offsetX
        }
        const target = element(event)
        const bounds = target.getBoundingClientRect()
        return (event as any).clientX - bounds.left  // 'any' to support older browsers
    }
    return 0
}

export function y(event: MouseEvent) {
    if (typeof event === 'object') {
        if ('offsetY' in event) {
            return event.offsetY
        }
        const target = element(event)
        const bounds = target.getBoundingClientRect()
        return (event as any).clientY - bounds.top  // 'any' to support older browsers
    }
    return 0
}