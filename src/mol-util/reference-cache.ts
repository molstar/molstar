/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export interface Reference<T> { value: T, usageCount: number }

export function createReference<T>(value: T, usageCount = 0) {
    return { value, usageCount }
}

export interface ReferenceItem<T> {
    free: () => void
    value: T
}

export function createReferenceItem<T>(ref: Reference<T>) {
    return {
        free: () => {
            ref.usageCount -= 1
        },
        value: ref.value
    }
}

export interface ReferenceCache<T, P, C> {
    get: (ctx: C, props: P) => ReferenceItem<T>
    clear: () => void
    count: number
    dispose: () => void
}

export function createReferenceCache<T, P, C>(hashFn: (props: P) => string, ctor: (ctx: C, props: P) => T, deleteFn: (v: T) => void): ReferenceCache<T, P, C> {
    const map: Map<string, Reference<T>> = new Map()

    return {
        get: (ctx: C, props: P) => {
            const id = hashFn(props)
            let ref = map.get(id)
            if (!ref) {
                ref = createReference<T>(ctor(ctx, props))
                map.set(id, ref)
            }
            ref.usageCount += 1
            return createReferenceItem(ref)
        },
        clear: () => {
            map.forEach((ref, id) => {
                if (ref.usageCount <= 0) {
                    if (ref.usageCount < 0) {
                        console.warn('Reference usageCount below zero.')
                    }
                    deleteFn(ref.value)
                }
            })
        },
        get count () {
            return map.size
        },
        dispose: () => {
            map.forEach(ref => deleteFn(ref.value))
        },
    }
}