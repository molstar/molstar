/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject } from 'rxjs';
import { MVSData } from '../../extensions/mvs/mvs-data';
import type { MVSStoriesViewerModel } from './elements/viewer';

export type MVSStoriesCommand =
    | { kind: 'load-mvs', format?: 'mvsj' | 'mvsx', url?: string, data?: MVSData | string | Uint8Array }


export class MVSStoriesContext {
    commands = new BehaviorSubject<MVSStoriesCommand | undefined>(undefined);
    state = {
        viewers: new BehaviorSubject<{ name?: string, model: MVSStoriesViewerModel }[]>([]),
        currentStoryData: new BehaviorSubject<string | Uint8Array | undefined>(undefined),
        isLoading: new BehaviorSubject(false),
    };

    dispatch(command: MVSStoriesCommand) {
        this.commands.next(command);
    }

    constructor(public name?: string) {
    }
}

export function getMVSStoriesContext(options?: { name?: string, container?: object }): MVSStoriesContext {
    const container: any = options?.container ?? window;
    container.componentContexts ??= {};
    const name = options?.name ?? '<default>';
    if (!container.componentContexts[name]) {
        container.componentContexts[name] = new MVSStoriesContext(options?.name);
    }
    return container.componentContexts[name];
}