

import { loadCheckpoint } from '../../mol-util/debug';
loadCheckpoint(`mol-task/execution/synchronous.ts::start`);
/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RuntimeContext } from './runtime-context';

export class SynchronousRuntimeContext implements RuntimeContext {
    shouldUpdate = false;
    isSynchronous = true;
    update(progress: string | Partial<RuntimeContext.ProgressUpdate>, dontNotify?: boolean): Promise<void> | void { }
}

export const SyncRuntimeContext = new SynchronousRuntimeContext();
loadCheckpoint(`mol-task/execution/synchronous.ts::end`);
