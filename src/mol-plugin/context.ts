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

export class PluginContext {
    private disposed = false;
    private ev = RxEventHelper.create();

    readonly state = new PluginState(this);
    readonly commands = new PluginCommand.Manager();

    readonly events = {
        state: {
            data: this.state.data.events,
            behavior: this.state.behavior.events
        }
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
        BuiltInPluginBehaviors.State.ApplyAction(this);
        BuiltInPluginBehaviors.State.RemoveObject(this);
        BuiltInPluginBehaviors.State.SetCurrentObject(this);
        BuiltInPluginBehaviors.State.Update(this);
    }

    async _test_initBehaviors() {
        const tree = this.state.behavior.tree.build()
            .toRoot().apply(PluginBehaviors.Representation.AddRepresentationToCanvas, { ref: PluginBehaviors.Representation.AddRepresentationToCanvas.id })
            .and().toRoot().apply(PluginBehaviors.Representation.HighlightLoci, { ref: PluginBehaviors.Representation.HighlightLoci.id })
            .and().toRoot().apply(PluginBehaviors.Representation.SelectLoci, { ref: PluginBehaviors.Representation.SelectLoci.id })
            .getTree();

        await this.runTask(this.state.behavior.update(tree));
    }

    _test_initDataActions() {
        this.state.data.actions
            .add(CreateStructureFromPDBe)
            .add(StateTransforms.Data.Download.toAction());
    }

    applyTransform(state: State, a: Transform.Ref, transformer: Transformer, params: any) {
        const tree = state.tree.build().to(a).apply(transformer, params);
        return PluginCommands.State.Update.dispatch(this, { state, tree });
    }

    updateTransform(state: State, a: Transform.Ref, params: any) {
        const tree = state.build().to(a).update(params);
        return PluginCommands.State.Update.dispatch(this, { state, tree });
    }

    _test_createState(id: string) {
        this.runTask(this.state.data.apply(CreateStructureFromPDBe, { id }));
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

        merge(this.events.state.data.object.replaced, this.events.state.behavior.object.replaced).subscribe(o => {
            if (o.oldObj && SO.isBehavior(o.oldObj)) o.oldObj.data.unregister();
            if (o.newObj && SO.isBehavior(o.newObj)) o.newObj.data.register();
        });
    }

    _test_centerView() {
        const sel = this.state.data.select(q => q.root.subtree().ofType(SO.Molecule.Structure.type));
        if (!sel.length) return;

        const center = (sel[0].obj! as SO.Molecule.Structure).data.boundary.sphere.center;
        this.canvas3d.camera.setState({ target: center });
        this.canvas3d.requestDraw(true);
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