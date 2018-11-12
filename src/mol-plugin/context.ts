/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Transformer, Transform, State } from 'mol-state';
import { Canvas3D } from 'mol-canvas3d/canvas3d';
import { StateTransforms } from './state/transforms';
import { PluginStateObject as SO } from './state/objects';
import { RxEventHelper } from 'mol-util/rx-event-helper';
import { PluginState } from './state';
import { PluginCommand, PluginCommands } from './command';
import { Task } from 'mol-task';
import { merge } from 'rxjs';
import { PluginBehaviors, BuiltInPluginBehaviors } from './behavior';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { Representation } from 'mol-repr';
import { CreateStructureFromPDBe } from './state/actions/basic';
import { LogEntry } from 'mol-util/log-entry';

export class PluginContext {
    private disposed = false;
    private ev = RxEventHelper.create();

    readonly state = new PluginState(this);
    readonly commands = new PluginCommand.Manager();

    readonly events = {
        state: {
            data: this.state.data.events,
            behavior: this.state.behavior.events,
            cameraSnapshots: this.state.cameraSnapshots.events,
            snapshots: this.state.snapshots.events,
        },
        log: this.ev<LogEntry>()
    };

    readonly behaviors = {
        state: {
            data: this.state.data.behaviors,
            behavior: this.state.behavior.behaviors
        },
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
            console.log('canvas3d created');
            return true;
        } catch (e) {
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

    async runTask<T>(task: Task<T>) {
        return await task.run(p => console.log(p.root.progress.message), 250);
    }

    dispose() {
        if (this.disposed) return;
        this.commands.dispose();
        this.canvas3d.dispose();
        this.ev.dispose();
        this.state.dispose();
        this.disposed = true;
    }

    private initBuiltInBehavior() {
        BuiltInPluginBehaviors.State.registerDefault(this);
        BuiltInPluginBehaviors.Representation.registerDefault(this);
        BuiltInPluginBehaviors.Camera.registerDefault(this);

        merge(this.state.data.events.log, this.state.behavior.events.log).subscribe(e => this.events.log.next(e));
    }

    async _test_initBehaviors() {
        const tree = this.state.behavior.tree.build()
            .toRoot().apply(PluginBehaviors.Representation.HighlightLoci, { ref: PluginBehaviors.Representation.HighlightLoci.id })
            .toRoot().apply(PluginBehaviors.Representation.SelectLoci, { ref: PluginBehaviors.Representation.SelectLoci.id })
            .getTree();

        await this.runTask(this.state.behavior.update(tree));
    }

    _test_initDataActions() {
        this.state.data.actions
            .add(CreateStructureFromPDBe)
            .add(StateTransforms.Data.Download)
            .add(StateTransforms.Model.CreateStructureAssembly)
            .add(StateTransforms.Model.CreateStructure)
            .add(StateTransforms.Model.CreateModelFromTrajectory)
            .add(StateTransforms.Visuals.CreateStructureRepresentation);
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
        merge(this.events.state.data.object.created, this.events.state.behavior.object.created).subscribe(o => {
            if (!SO.isBehavior(o.obj)) return;
            console.log('registering behavior', o.obj.label);
            o.obj.data.register();
        });

        merge(this.events.state.data.object.removed, this.events.state.behavior.object.removed).subscribe(o => {
            if (!SO.isBehavior(o.obj)) return;
            o.obj.data.unregister();
        });

        merge(this.events.state.data.object.updated, this.events.state.behavior.object.updated).subscribe(o => {
            if (o.action === 'recreate') {
                if (o.oldObj && SO.isBehavior(o.oldObj)) o.oldObj.data.unregister();
                if (o.obj && SO.isBehavior(o.obj)) o.obj.data.register();
            }
        });
    }

    constructor() {
        this.initEvents();
        this.initBuiltInBehavior();

        this._test_initBehaviors();
        this._test_initDataActions();
    }

    // logger = ;
    // settings = ;
}