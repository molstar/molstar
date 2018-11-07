/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTree, StateSelection, Transformer, Transform } from 'mol-state';
import Canvas3D from 'mol-canvas3d/canvas3d';
import { StateTransforms } from './state/transforms';
import { PluginStateObjects as SO } from './state/objects';
import { RxEventHelper } from 'mol-util/rx-event-helper';
import { PluginState } from './state';
import { MolScriptBuilder } from 'mol-script/language/builder';
import { PluginCommand, PluginCommands } from './command';
import { Task } from 'mol-task';
import { merge } from 'rxjs';
import { PluginBehaviors } from './behavior';

export class PluginContext {
    private disposed = false;
    private ev = RxEventHelper.create();

    readonly state = new PluginState(this);
    readonly commands = new PluginCommand.Manager();

    readonly events = {
        state: {
            data: this.state.data.context.events,
            behavior: this.state.behavior.context.events
        }
    };

    readonly behaviors = {
        state: {
            data: this.state.data.context.behaviors,
            behavior: this.state.behavior.context.behaviors
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
        const req = await fetch(url);
        return type === 'string' ? await req.text() : new Uint8Array(await req.arrayBuffer());
    }

    async runTask<T>(task: Task<T>) {
        return await task.run(p => console.log(p), 250);
    }

    dispose() {
        if (this.disposed) return;
        this.commands.dispose();
        this.canvas3d.dispose();
        this.ev.dispose();
        this.state.dispose();
        this.disposed = true;
    }

    async _test_initBehaviours() {
        const tree = StateTree.build(this.state.behavior.tree)
            .toRoot().apply(PluginBehaviors.Data.SetCurrentObject)
            .and().toRoot().apply(PluginBehaviors.Data.Update)
            .and().toRoot().apply(PluginBehaviors.Data.RemoveObject)
            .and().toRoot().apply(PluginBehaviors.Representation.AddRepresentationToCanvas)
            .getTree();

        await this.state.updateBehaviour(tree);
    }

    _test_applyTransform(a: Transform.Ref, transformer: Transformer, params: any) {
        const tree = StateTree.build(this.state.data.tree).to(a).apply(transformer, params).getTree();
        PluginCommands.Data.Update.dispatch(this, { tree });
    }

    _test_createState(url: string) {
        const b = StateTree.build(this.state.data.tree);

        const query = MolScriptBuilder.struct.generator.atomGroups({
            // 'atom-test': MolScriptBuilder.core.rel.eq([
            //     MolScriptBuilder.struct.atomProperty.macromolecular.label_comp_id(),
            //     MolScriptBuilder.es('C')
            // ]),
            'residue-test': MolScriptBuilder.core.rel.eq([
                MolScriptBuilder.struct.atomProperty.macromolecular.label_comp_id(),
                'ALA'
            ])
        });

        const newTree = b.toRoot()
            .apply(StateTransforms.Data.Download, { url })
            .apply(StateTransforms.Data.ParseCif)
            .apply(StateTransforms.Model.ParseModelsFromMmCif, {}, { ref: 'models' })
            .apply(StateTransforms.Model.CreateStructureFromModel, { modelIndex: 0 }, { ref: 'structure' })
            .apply(StateTransforms.Model.CreateStructureAssembly)
            .apply(StateTransforms.Model.CreateStructureSelection, { query, label: 'ALA residues' })
            .apply(StateTransforms.Visuals.CreateStructureRepresentation)
            .getTree();

        this.state.updateData(newTree);
    }

    private initEvents() {
        merge(this.events.state.data.object.created, this.events.state.behavior.object.created).subscribe(o => {
            console.log('creating', o.obj.type);
            if (!SO.Behavior.is(o.obj)) return;
            o.obj.data.register();
        });

        merge(this.events.state.data.object.removed, this.events.state.behavior.object.removed).subscribe(o => {
            if (!SO.Behavior.is(o.obj)) return;
            o.obj.data.unregister();
        });

        merge(this.events.state.data.object.replaced, this.events.state.behavior.object.replaced).subscribe(o => {
            if (o.oldObj && SO.Behavior.is(o.oldObj)) o.oldObj.data.unregister();
            if (o.newObj && SO.Behavior.is(o.newObj)) o.newObj.data.register();
        });
    }

    _test_centerView() {
        const sel = StateSelection.select(StateSelection.root().subtree().ofType(SO.Structure.type), this.state.data);
        if (!sel.length) return;

        const center = (sel[0].obj! as SO.Structure).data.boundary.sphere.center;
        console.log({ sel, center, rc: this.canvas3d.reprCount });
        this.canvas3d.center(center);
        this.canvas3d.requestDraw(true);
    }

    _test_nextModel() {
        const models = StateSelection.select('models', this.state.data)[0].obj as SO.Models;
        const idx = (this.state.data.tree.getValue('structure')!.params as Transformer.Params<typeof StateTransforms.Model.CreateStructureFromModel>).modelIndex;
        const newTree = StateTree.updateParams(this.state.data.tree, 'structure', { modelIndex: (idx + 1) % models.data.length });
        return this.state.updateData(newTree);
        // this.viewer.requestDraw(true);
    }

    _test_playModels() {
        const update = async () => {
            await this._test_nextModel();
            setTimeout(update, 1000 / 15);
        }
        update();
    }

    constructor() {
        this.initEvents();

        this._test_initBehaviours();
    }

    // logger = ;
    // settings = ;
}