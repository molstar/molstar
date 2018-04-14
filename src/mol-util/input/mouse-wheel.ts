/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/mikolalysenko/mouse-wheel,
 * copyright (c) 2015 Mikola Lysenko. MIT License
 */

import { Subject } from 'rxjs';
import toPixels from '../to-pixels'

interface MouseWheel {
    noScroll: boolean
    wheel: Subject<{ dx: number, dy: number, dz: number, event: WheelEvent }>
    dispose: () => void
}

namespace MouseWheel {
    export function create(element: Element, noScroll = true): MouseWheel {
        const lineHeight = toPixels('ex', element)
        let disposed = false
        const wheel = new Subject<{ dx: number, dy: number, dz: number, event: WheelEvent }>()

        element.addEventListener('wheel', listener)

        return {
            get noScroll () { return noScroll },
            set noScroll (value: boolean) { noScroll = value },

            wheel,
            dispose
        }

        function listener(event: MouseWheelEvent) {
            if (noScroll) {
                event.preventDefault()
            }
            const mode = event.deltaMode
            let dx = event.deltaX || 0
            let dy = event.deltaY || 0
            let dz = event.deltaZ || 0
            let scale = 1
            switch (mode) {
                case 1: scale = lineHeight; break
                case 2: scale = window.innerHeight; break
            }
            dx *= scale
            dy *= scale
            dz *= scale
            if (dx || dy || dz) {
                wheel.next({ dx, dy, dz, event })
            }
        }

        function dispose() {
            if (disposed) return
            disposed = true
            element.removeEventListener('wheel', listener)
            wheel.unsubscribe()
        }
    }
}

export default MouseWheel