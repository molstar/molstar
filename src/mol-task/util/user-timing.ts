/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '../task';
import { isProductionMode } from '../../mol-util/debug';

const hasPerformance = (typeof performance !== 'undefined') && performance.mark && performance.measure;
const timingEnabled = hasPerformance && !isProductionMode;

export namespace UserTiming {
    function startMarkName(task: Task<any>) { return `startTask${task.id}`; }
    function endMarkName(task: Task<any>) { return `endTask${task.id}`; }
    export function markStart(task: Task<any>) {
        if (timingEnabled) performance.mark(startMarkName(task));
    }
    export function markEnd(task: Task<any>) {
        if (timingEnabled) performance.mark(endMarkName(task));
    }
    export function measure(task: Task<any>) {
        if (timingEnabled) performance.measure(`✳️ ${task.name}`, startMarkName(task), endMarkName(task));
    }
}