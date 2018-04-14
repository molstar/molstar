/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable } from './renderable';

export default function createStats (renderables: Renderable[]) {
    const prevGpuTimes: number[] = []
    for (let i = 0; i < renderables.length; i++) {
        prevGpuTimes[i] = 0
    }

    let frameTimeCount = 0
    let totalTime = 1.1
    let N = 50

    const totalFrameTime: number[] = []
    const avgFrameTime: number[] = []
    for (let i = 0; i < renderables.length; ++i) {
        totalFrameTime[i] = 0.0
        avgFrameTime[i] = 0.0
    }

    return {
        add: (renderable: Renderable) => {
            renderables.push(renderable)
            prevGpuTimes.push(0)
            totalFrameTime.push(0)
            avgFrameTime.push(0)
        },
        update: (deltaTime: number) => {
            totalTime += deltaTime
            if (totalTime > 1.0) {
                totalTime = 0

                // for (let i = 0; i < renderables.length; i++) {
                //     const renderable = renderables[i]
                //     const str = `${renderable.name}: ${Math.round(100.0 * avgFrameTime[i]) / 100.0}ms`
                //     console.log(str)
                // }

                const sumFrameTime = avgFrameTime.reduce((x: number, y: number) => x + y, 0)
                const str = `${Math.round(100.0 * sumFrameTime) / 100.0}ms`
                console.log(str)
            }

            frameTimeCount++

            for (let i = 0; i < renderables.length; i++) {
                const renderable = renderables[i]
                const frameTime = renderable.stats.gpuTime - prevGpuTimes[i]
                totalFrameTime[i] += frameTime

                if (frameTimeCount === N) {
                    avgFrameTime[i] = totalFrameTime[i] / N
                    totalFrameTime[i] = 0.0
                }

                prevGpuTimes[i] = renderable.stats.gpuTime
            }

            if (frameTimeCount === N) frameTimeCount = 0
        }
    }
}