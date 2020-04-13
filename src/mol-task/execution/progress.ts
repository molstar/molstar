/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '../task';

interface Progress {
    root: Progress.Node,
    canAbort: boolean,
    requestAbort: (reason?: string) => void
}

namespace Progress {
    export interface Node {
        readonly progress: Task.Progress,
        readonly children: ReadonlyArray<Node>
    }

    export interface Observer { (progress: Progress): void }

    function _format(root: Progress.Node, prefix = ''): string {
        const p = root.progress;
        if (!root.children.length) {
            if (p.isIndeterminate) return `${prefix}${p.taskName}: ${p.message}`;
            return `${prefix}${p.taskName}: [${p.current}/${p.max}] ${p.message}`;
        }

        const newPrefix = prefix + '  |_ ';
        const subTree = root.children.map(c => _format(c, newPrefix));
        if (p.isIndeterminate) return `${prefix}${p.taskName}: ${p.message}\n${subTree.join('\n')}`;
        return `${prefix}${p.taskName}: [${p.current}/${p.max}] ${p.message}\n${subTree.join('\n')}`;
    }

    export function format(p: Progress) { return _format(p.root); }
}

export { Progress };