/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getCellBoundingSphere } from '../../mol-plugin-state/manager/focus-camera/focus-object';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateObjectCell } from '../../mol-state';
import { PluginContext } from '../../mol-plugin/context';
import { Script } from '../../mol-script/script';
import { QueryContext, QueryFn, StructureElement, StructureSelection } from '../../mol-model/structure';
import { BehaviorSubject } from 'rxjs';

export type MarkdownExtensionEvent = 'click' | 'mouse-enter' | 'mouse-leave';

export interface MarkdownExtension {
    name: string;
    execute?: (options: {
        event: MarkdownExtensionEvent,
        args: Record<string, string>,
        manager: MarkdownExtensionManager
    }) => void;
    reactRenderFn?: (options: {
        args: Record<string, string>,
        manager: MarkdownExtensionManager
    }) => any;
}

export const BuiltInMarkdownExtension: MarkdownExtension[] = [
    {
        name: 'center-camera',
        execute: ({ event, args, manager }) => {
            if (event !== 'click') return;
            if ('center-camera' in args) {
                manager.plugin.managers.camera.reset();
            }
        }
    },
    {
        name: 'apply-snapshot',
        execute: ({ event, args, manager }) => {
            if (event !== 'click') return;
            const key = args['apply-snapshot'];
            if (!key) return;
            manager.plugin.managers.snapshot.applyKey(key);
        }
    },
    {
        name: 'focus-refs',
        execute: ({ event, args, manager }) => {
            if (event !== 'click') return;
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
        execute: ({ event, args, manager }) => {
            const refs = parseArray(args['highlight-refs']);
            if (!refs?.length) return;

            if (event === 'mouse-leave' && refs.length) {
                manager.plugin.managers.interactivity.lociHighlights.clearHighlights();
                return;
            } else if (event === 'mouse-enter') {
                const cells = manager.findCells(refs);
                for (const cell of findRepresentations(manager.plugin, cells)) {
                    if (!cell.obj?.data) continue;
                    const { repr } = cell.obj.data;
                    for (const loci of repr.getAllLoci()) {
                        manager.plugin.managers.interactivity.lociHighlights.highlight({ loci, repr }, false);
                    }
                }
            }
        }
    },
    {
        name: 'query',
        execute: ({ event, args, manager }) => {
            const expression = args['query'];
            if (!expression?.length) return;

            // supported languages: mol-script, pymol, vmd, jmol
            const language = args['lang'] || 'mol-script';
            // supported actions: highlight, focus
            const action = parseArray(args['action'] || 'highlight');
            const focusRadius = parseFloat(args['focus-radius'] || '3');

            if (event === 'mouse-leave') {
                if (action.includes('highlight')) {
                    manager.plugin.managers.interactivity.lociHighlights.clearHighlights();
                }
                return;
            }

            let query: QueryFn<StructureSelection>;
            try {
                query = Script.toQuery({
                    language: language as Script.Language,
                    expression
                });
            } catch (e) {
                console.warn(`Failed to parse query '${expression}' (${language})`, e);
                return;
            }

            const structures = manager.plugin.state.data.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Structure));

            if (event === 'mouse-enter') {
                if (!action.includes('focus')) {
                    return;
                }
                manager.plugin.managers.interactivity.lociHighlights.clearHighlights();
                for (const structure of structures) {
                    if (!structure.obj?.data) continue;
                    const selection = query(new QueryContext(structure.obj.data));
                    const loci = StructureSelection.toLociWithSourceUnits(selection);
                    manager.plugin.managers.interactivity.lociHighlights.highlight({
                        loci,
                    }, false);
                }
            }

            if (event === 'click') {
                if (!action.includes('focus')) {
                    return;
                }
                const spheres = structures.map(s => {
                    if (!s.obj?.data) return undefined;
                    const selection = query(new QueryContext(s.obj.data));
                    if (StructureSelection.isEmpty(selection)) return;

                    const loci = StructureSelection.toLociWithSourceUnits(selection);
                    return StructureElement.Loci.getBoundary(loci).sphere;
                }).filter(s => !!s);

                if (spheres.length) {
                    manager.plugin.managers.camera.focusSpheres(spheres, s => s, { extraRadius: focusRadius });
                }
            }
        },
    },
    {
        name: 'play-audio',
        execute: ({ event, args, manager }) => {
            if (event !== 'click') return;

            const src = args['play-audio'];
            if (!src?.length) return;
            manager.audio.play(src);
        }
    },
    {
        name: 'toggle-audio',
        execute: ({ event, args, manager }) => {
            if (event !== 'click' || !('toggle-audio' in args)) return;

            const src = args['toggle-audio'];
            manager.audio.play(src, { toggle: true });
        }
    },
    {
        name: 'pause-audio',
        execute: ({ event, args, manager }) => {
            if (event !== 'click' || !('pause-audio' in args)) return;
            manager.audio.pause();
        }
    },
    {
        name: 'stop-audio',
        execute: ({ event, args, manager }) => {
            if (event !== 'click' || !('stop-audio' in args)) return;
            manager.audio.stop();
        }
    },
];

export class MarkdownExtensionManager {
    state = {
        audioPlayer: new BehaviorSubject<HTMLAudioElement | null>(null),
    };

    private extension: MarkdownExtension[] = [];
    private refResolvers: Record<string, (plugin: PluginContext, refs: string[]) => StateObjectCell[]> = {
        default: (plugin: PluginContext, refs: string[]) => refs
            .map(ref => plugin.state.data.cells.get(ref))
            .filter(c => !!c),
    };
    private uriResolvers: Record<string, (plugin: PluginContext, uri: string) => Promise<string> | string | undefined> = {};
    private argsParsers: [name: string, priority: number, parser: (input: string | undefined) => Record<string, string> | undefined][] = [
        ['default', 100, defaultParseMarkdownCommandArgs],
    ];

    /**
     * Default parser has priority 100, parsers with higher priority
     * will be called first.
     */
    registerArgsParser(name: string, priority: number, parser: (input: string | undefined) => Record<string, string> | undefined) {
        this.removeArgsParser(name);
        this.argsParsers.push([name, priority, parser]);
        this.argsParsers.sort((a, b) => b[1] - a[1]); // Sort by priority, higher first
    }

    removeArgsParser(name: string) {
        const idx = this.argsParsers.findIndex(p => p[0] === name);
        if (idx >= 0) {
            this.argsParsers.splice(idx, 1);
        }
    }

    parseArgs(input: string | undefined): Record<string, string> | undefined {
        for (const [,, parser] of this.argsParsers) {
            const ret = parser(input);
            if (ret) return ret;
        }
        return undefined;
    }

    registerRefResolver(name: string, resolver: (plugin: PluginContext, refs: string[]) => StateObjectCell[]) {
        this.refResolvers[name] = resolver;
    }

    removeRefResolver(name: string) {
        delete this.refResolvers[name];
    }

    registerUriResolver(name: string, resolver: (plugin: PluginContext, uri: string) => Promise<string> | string | undefined) {
        this.uriResolvers[name] = resolver;
    }

    removeUriResolver(name: string) {
        delete this.uriResolvers[name];
    }

    registerExtension(command: MarkdownExtension) {
        const existing = this.extension.findIndex(c => c.name === command.name);
        if (existing >= 0) {
            this.extension[existing] = command;
        } else {
            this.extension.push(command);
        }
    }

    removeExtension(name: string) {
        const idx = this.extension.findIndex(c => c.name === name);
        if (idx >= 0) {
            this.extension.splice(idx, 1);
        }
    }

    private _tryRender(ext: MarkdownExtension, options: { args: Record<string, string>, manager: MarkdownExtensionManager }) {
        try {
            return ext.reactRenderFn?.(options);
        } catch (e) {
            console.error(`Failed to render markdown extension '${ext.name}'`, e);
            return null;
        }
    }

    /**
     * Render a markdown extension with the given arguments.
     * Default renderers are defined separately because we
     * don't want to include `react` outside of mol-plugin-ui.
     */
    tryRender(args: Record<string, string>, defaultRenderers: MarkdownExtension[]): any {
        const options = { args, manager: this };
        for (const ext of this.extension) {
            const ret = this._tryRender(ext, options);
            if (ret) {
                return ret;
            }
        }
        for (const ext of defaultRenderers) {
            const ret = this._tryRender(ext, options);
            if (ret) {
                return ret;
            }
        }
        return null;
    }

    tryExecute(event: MarkdownExtensionEvent, args: Record<string, string>) {
        const options = { event, args, manager: this };
        for (const ext of this.extension) {
            try {
                ext.execute?.(options);
            } catch (e) {
                console.error(`Failed to execute markdown extension '${ext.name}'`, e);
            }
        }
    }

    tryResolveUri(uri: string): Promise<string> | string | undefined {
        for (const resolver of Object.values(this.uriResolvers)) {
            const resolved = resolver(this.plugin, uri);
            if (resolved) {
                return resolved;
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

    private resolveAudioPlayer() {
        if (this.state.audioPlayer.value) {
            return this.state.audioPlayer.value;
        }

        const audio = document.createElement('audio');
        audio.controls = true;
        audio.preload = 'auto';
        audio.style.width = '100%';
        audio.style.height = '32px';
        this.state.audioPlayer.next(audio);
        return audio;
    }

    get audioPlayer() {
        return this.state.audioPlayer.value;
    }

    audio = {
        play: async (src: string, options?: { toggle?: boolean }) => {
            try {
                const audio = this.resolveAudioPlayer();

                let newSource = false;
                if (src?.trim()) {
                    const resolved = this.tryResolveUri(src);
                    let uri: string = src;
                    if (typeof (resolved as Promise<string>)?.then === 'function') {
                        uri = (await resolved) as string;
                    } else if (resolved) {
                        uri = resolved as string;
                    }
                    newSource = audio.src !== uri;
                    if (newSource) {
                        audio.src = uri;
                        audio.load();
                    }
                }

                if (!newSource && options?.toggle) {
                    if (audio.paused) {
                        await audio.play();
                    } else {
                        audio.pause();
                    }
                } else {
                    audio.currentTime = 0;
                    await audio.play();
                }
            } catch (e) {
                console.error('Failed to play audio', e);
            }
        },
        pause: () => {
            this.audioPlayer?.pause();
        },
        stop: () => {
            if (!this.audioPlayer) return;
            this.audioPlayer.pause();
            this.audioPlayer.currentTime = 0;
        },
        dispose: () => {
            if (this.audioPlayer) {
                this.audioPlayer.pause();
                this.audioPlayer.currentTime = 0;
                this.state.audioPlayer.next(null);
            }
        }
    };

    constructor(public plugin: PluginContext) {
        for (const command of BuiltInMarkdownExtension) {
            this.registerExtension(command);
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

export function defaultParseMarkdownCommandArgs(input: string | undefined): Record<string, string> | undefined {
    if (!input?.startsWith('!')) return undefined;
    const entries = decodeURIComponent(input.substring(1))
        .split('&')
        .map(p => p.trim())
        .filter(p => p.length > 0)
        .map(p => p.split('=', 2).map(s => s.trim()));
    if (entries.length === 0) return undefined;
    return Object.fromEntries(entries);
}