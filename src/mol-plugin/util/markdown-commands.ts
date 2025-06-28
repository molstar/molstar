/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getCellBoundingSphere } from '../../mol-plugin-state/manager/focus-camera/focus-object';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateObjectCell } from '../../mol-state';
import { PluginContext } from '../context';

export interface MarkdownCommand {
    name: string;
    event: 'click' | 'mouse-enter' | 'mouse-leave';
    execute: (args: Record<string, string>, manager: MarkdownCommandManager) => void;
}

export interface MarkdownRenderer {
    name: string;
    reactRenderFn: (args: Record<string, string>, manager: MarkdownCommandManager) => any;
}

export const DefaultMarkdownCommands: MarkdownCommand[] = [
    {
        name: 'center-camera',
        event: 'click',
        execute: (args, manager) => {
            if ('center-camera' in args) {
                manager.plugin.managers.camera.reset();
            }
        }
    },
    {
        name: 'apply-snapshot',
        event: 'click',
        execute: (args, manager) => {
            const key = args['apply-snapshot'];
            if (!key) return;
            manager.plugin.managers.snapshot.applyKey(key);
        }
    },
    {
        name: 'focus-refs',
        event: 'click',
        execute: (args, manager) => {
            const refs = parseArray(args['focus-refs']);
            if (!refs?.length) return;

            const cells = manager.findCells(refs);
            if (!cells.length) return;

            const reprs = findRepresentations(manager.plugin, cells);
            if (!reprs.length) return;

            const spheres = reprs.map(c => getCellBoundingSphere(manager.plugin, c.transform.ref)).filter(s => !!s);
            if (!spheres.length) return;
            manager.plugin.managers.camera.focusSpheres(spheres, s => s, { extraRadius: 3 });
        }
    },
    {
        name: 'highlight-refs',
        event: 'mouse-enter',
        execute: (args, manager) => {
            const refs = parseArray(args['highlight-refs']);
            if (!refs?.length) return;

            const cells = manager.findCells(refs);
            for (const cell of findRepresentations(manager.plugin, cells)) {
                if (!cell.obj?.data) continue;
                const { repr } = cell.obj.data;
                for (const loci of repr.getAllLoci()) {
                    manager.plugin.managers.interactivity.lociHighlights.highlight({ loci, repr }, false);
                }
            }
        }
    },
    {
        name: 'highlight-refs',
        event: 'mouse-leave',
        execute: (args, manager) => {
            if ('highlight-refs' in args) {
                manager.plugin.managers.interactivity.lociHighlights.clearHighlights();
            }
        }
    }
];

export class MarkdownCommandManager {
    commands: MarkdownCommand[] = [];
    renderers: MarkdownRenderer[] = [];

    refResolvers: Record<string, (plugin: PluginContext, refs: string[]) => StateObjectCell[]> = {
        default: (plugin: PluginContext, refs: string[]) => refs
            .map(ref => plugin.state.data.cells.get(ref))
            .filter(c => !!c),
    };

    registerCommand(command: MarkdownCommand) {
        const existing = this.commands.findIndex(c => c.name === command.name && c.event === command.event);
        if (existing >= 0) {
            this.commands[existing] = command;
        } else {
            this.commands.push(command);
        }
    }

    removeCommand(name: string, event: MarkdownCommand['event']) {
        const idx = this.commands.findIndex(c => c.name === name && c.event === event);
        if (idx >= 0) {
            this.commands.splice(idx, 1);
        }
    }

    registerRenderer(renderer: MarkdownRenderer) {
        const existing = this.renderers.findIndex(r => r.name === renderer.name);
        if (existing >= 0) {
            this.renderers[existing] = renderer;
        } else {
            this.renderers.push(renderer);
        }
    }

    removeRenderer(name: string) {
        const idx = this.renderers.findIndex(r => r.name === name);
        if (idx >= 0) {
            this.renderers.splice(idx, 1);
        }
    }

    render(args: Record<string, string>, defaultRenderers: MarkdownRenderer[]): any {
        for (const renderer of this.renderers) {
            const ret = renderer.reactRenderFn(args, this);
            if (ret) {
                return ret;
            }
        }
        for (const renderer of defaultRenderers) {
            const ret = renderer.reactRenderFn(args, this);
            if (ret) {
                return ret;
            }
        }
        return null;
    }

    execute(event: MarkdownCommand['event'], args: Record<string, string>) {
        for (const command of this.commands) {
            if (command.event === event) {
                command.execute(args, this);
            }
        }
    }

    findCells(refs: string[]): StateObjectCell[] {
        const added = new Set<string>();
        const ret: StateObjectCell[] = [];
        for (const resolver of Object.values(this.refResolvers)) {
            for (const cell of resolver(this.plugin, refs)) {
                if (cell && !added.has(cell.transform.ref)) {
                    added.add(cell.transform.ref);
                    ret.push(cell);
                }
            }
        }
        return ret;
    }

    constructor(public plugin: PluginContext) {
        for (const command of DefaultMarkdownCommands) {
            this.registerCommand(command);
        }
    }
}

function parseArray(input?: string): string[] {
    return input?.split(',').map(s => s.trim()).filter(s => s.length > 0) ?? [];
}

function findRepresentations(plugin: PluginContext, cells: StateObjectCell[]): StateObjectCell[] {
    if (!cells.length) return [];
    return plugin.state.data.selectQ(q =>
        q.byValue(...cells).subtree().filter(c => PluginStateObject.isRepresentation3D(c.obj))
    );
}

export function parseMarkdownCommandArgs(input: string): Record<string, string> {
    return Object.fromEntries(decodeURIComponent(input)
        .split('&')
        .map(p => p.trim())
        .filter(p => p.length > 0)
        .map(p => p.split('=', 2).map(s => s.trim())));
}