/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '../task'

const hasPerformance = typeof performance !== 'undefined'
/**
 * on node `process.env.NODE_ENV` is available, in webpack build it is automatically set
 * by the DefinePlugin to the webpack `mode` value
 */
const productionMode = process.env.NODE_ENV === 'production'
const timingEnabled = hasPerformance && !productionMode

export namespace UserTiming {
    function startMarkName(task: Task<any>) { return `startTask${task.id}` }
    function endMarkName(task: Task<any>) { return `endTask${task.id}` }
    export function markStart(task: Task<any>) {
        if (timingEnabled) performance.mark(startMarkName(task))
    }
    export function markEnd(task: Task<any>) {
        if (timingEnabled) performance.mark(endMarkName(task))
    }
    export function measure(task: Task<any>) {
        if (timingEnabled) performance.measure(`✳️ ${task.name}`, startMarkName(task), endMarkName(task))
    }
}