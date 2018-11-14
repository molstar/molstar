/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Canvas3D } from 'mol-canvas3d/canvas3d';
import { EmptyLoci, Loci } from 'mol-model/loci';
import { Representation } from 'mol-repr/representation';
import { StructureRepresentationRegistry } from 'mol-repr/structure/registry';
import { State, Transform, Transformer } from 'mol-state';
import { Task } from 'mol-task';
import { ColorTheme } from 'mol-theme/color';
import { SizeTheme } from 'mol-theme/size';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { LogEntry } from 'mol-util/log-entry';
import { RxEventHelper } from 'mol-util/rx-event-helper';
import { merge } from 'rxjs';
import { BuiltInPluginBehaviors } from './behavior';
import { PluginCommand, PluginCommands } from './command';
import { PluginSpec } from './spec';
import { PluginState } from './state';
import { PluginStateObject as SO } from './state/objects';
import { TaskManager } from './util/task-manager';

export class PluginContext {
    private disposed = false;
    private ev = RxEventHelper.create();
    private tasks = new TaskManager();

    readonly state = new PluginState(this);
    readonly commands = new PluginCommand.Manager();

    readonly events = {
        state: {
            cell: {
                stateUpdated: merge(this.state.dataState.events.cell.stateUpdated, this.state.behaviorState.events.cell.stateUpdated),
                created: merge(this.state.dataState.events.cell.created, this.state.behaviorState.events.cell.created),
                removed: merge(this.state.dataState.events.cell.removed, this.state.behaviorState.events.cell.removed),
            },
            object: {
                created: merge(this.state.dataState.events.object.created, this.state.behaviorState.events.object.created),
                removed: merge(this.state.dataState.events.object.removed, this.state.behaviorState.events.object.removed),
                updated: merge(this.state.dataState.events.object.updated, this.state.behaviorState.events.object.updated)
            },
            // data: this.state.dataState.events,
            // behavior: this.state.behaviorState.events,
            cameraSnapshots: this.state.cameraSnapshots.events,
            snapshots: this.state.snapshots.events,
        },
        log: this.ev<LogEntry>(),
        task: this.tasks.events
    };

    readonly structureReprensentation = {
        registry: new StructureRepresentationRegistry(),
        themeCtx: { colorThemeRegistry: new ColorTheme.Registry(), sizeThemeRegistry: new SizeTheme.Registry() } as ThemeRegistryContext
    }

    readonly behaviors = {
        // state: {
        //     data: this.state.dataState.behaviors,
        //     behavior: this.state.behaviorState.behaviors
        // },
        canvas: {
            highlightLoci: this.ev.behavior<{ loci: Loci, repr?: Representation.Any }>({ loci: EmptyLoci }),
            selectLoci: this.ev.behavior<{ loci: Loci, repr?: Representation.Any }>({ loci: EmptyLoci }),
        },
        command: this.commands.behaviour
    };

    readonly canvas3d: Canvas3D;


    initViewer(canvas: HTMLCanvasElement, container: HTMLDivElement) {
        try {
            (this.canvas3d as Canvas3D) = Canvas3D.create(canvas, container);
            this.canvas3d.animate();
            return true;
        } catch (e) {
            this.log(LogEntry.error('' + e));
            console.error(e);
            return false;
        }
    }

    log(e: LogEntry) {
        this.events.log.next(e);
    }

    /**
     * This should be used in all transform related request so that it could be "spoofed" to allow
     * "static" access to resources.
     */
    async fetch(url: string, type: 'string' | 'binary' = 'string'): Promise<string | Uint8Array> {
        const req = await fetch(url, { referrerPolicy: 'origin-when-cross-origin' });
        return type === 'string' ? await req.text() : new Uint8Array(await req.arrayBuffer());
    }

    runTask<T>(task: Task<T>) {
        return this.tasks.run(task);
    }

    dispose() {
        if (this.disposed) return;
        this.commands.dispose();
        this.canvas3d.dispose();
        this.ev.dispose();
        this.state.dispose();
        this.tasks.dispose();
        this.disposed = true;
    }

    private initBuiltInBehavior() {
        BuiltInPluginBehaviors.State.registerDefault(this);
        BuiltInPluginBehaviors.Representation.registerDefault(this);
        BuiltInPluginBehaviors.Camera.registerDefault(this);

        merge(this.state.dataState.events.log, this.state.behaviorState.events.log).subscribe(e => this.events.log.next(e));
    }

    async initBehaviors() {
        const tree = this.state.behaviorState.tree.build();

        for (const b of this.spec.behaviors) {
            tree.toRoot().apply(b.transformer, b.defaultParams || { }, { ref: b.transformer.id });
        }

        await this.runTask(this.state.behaviorState.update(tree));
    }

    initDataActions() {
        for (const a of this.spec.actions) {
            this.state.dataState.actions.add(a.action);
        }
    }

    applyTransform(state: State, a: Transform.Ref, transformer: Transformer, params: any) {
        const tree = state.tree.build().to(a).apply(transformer, params);
        return PluginCommands.State.Update.dispatch(this, { state, tree });
    }

    updateTransform(state: State, a: Transform.Ref, params: any) {
        const tree = state.build().to(a).update(params);
        return PluginCommands.State.Update.dispatch(this, { state, tree });
    }

    private initEvents() {
        this.events.state.object.created.subscribe(o => {
            if (!SO.isBehavior(o.obj)) return;
            o.obj.data.register();
        });

        this.events.state.object.removed.subscribe(o => {
            if (!SO.isBehavior(o.obj)) return;
            o.obj.data.unregister();
        });

        this.events.state.object.updated.subscribe(o => {
            if (o.action === 'recreate') {
                if (o.oldObj && SO.isBehavior(o.oldObj)) o.oldObj.data.unregister();
                if (o.obj && SO.isBehavior(o.obj)) o.obj.data.register();
            }
        });
    }

    constructor(public spec: PluginSpec) {
        this.initEvents();
        this.initBuiltInBehavior();

        this.initBehaviors();
        this.initDataActions();
    }

    // logger = ;
    // settings = ;
}