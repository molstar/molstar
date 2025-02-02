/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject } from 'rxjs';
import { MVSData } from '../../extensions/mvs/mvs-data';
import type { MolComponentViewerModel } from './elements/viewer';

export type MolComponentCommand =
    | { kind: 'load-mvs', format?: 'mvsj' | 'mvsx', url?: string, data?: MVSData }


export class MolComponentContext {
    commands = new BehaviorSubject<MolComponentCommand | undefined>(undefined);
    behavior = {
        viewers: new BehaviorSubject<{ name?: string, model: MolComponentViewerModel }[]>([]),
    };

    dispatch(command: MolComponentCommand) {
        this.commands.next(command);
    }

    constructor(public name?: string) {
    }
}

export function getMolComponentContext(options?: { name?: string, container?: object }) {
    const container: any = options?.container ?? window;
    container.componentContexts ??= {};
    const name = options?.name ?? '<default>';
    if (!container.componentContexts[name]) {
        container.componentContexts[name] = new MolComponentContext(options?.name);
    }
    return container.componentContexts[name];
}