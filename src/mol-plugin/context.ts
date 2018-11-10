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
import { MolScriptBuilder } from 'mol-script/language/builder';
import { PluginCommand, PluginCommands } from './command';
import { Task } from 'mol-task';
import { merge } from 'rxjs';
import { PluginBehaviors } from './behavior';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { Representation } from 'mol-repr';

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
        const tree = this.state.behavior.tree.build()
            .toRoot().apply(PluginBehaviors.Data.SetCurrentObject, { ref: PluginBehaviors.Data.SetCurrentObject.id })
            .and().toRoot().apply(PluginBehaviors.Data.Update, { ref: PluginBehaviors.Data.Update.id })
            .and().toRoot().apply(PluginBehaviors.Data.RemoveObject, { ref: PluginBehaviors.Data.RemoveObject.id })
            .and().toRoot().apply(PluginBehaviors.Representation.AddRepresentationToCanvas, { ref: PluginBehaviors.Representation.AddRepresentationToCanvas.id })
            .and().toRoot().apply(PluginBehaviors.Representation.HighlightLoci, { ref: PluginBehaviors.Representation.HighlightLoci.id })
            .and().toRoot().apply(PluginBehaviors.Representation.SelectLoci, { ref: PluginBehaviors.Representation.SelectLoci.id })
            .getTree();

        await this.runTask(this.state.behavior.update(tree));
    }

    applyTransform(state: State, a: Transform.Ref, transformer: Transformer, params: any) {
        const tree = state.tree.build().to(a).apply(transformer, params);
        return PluginCommands.State.Update.dispatch(this, { state, tree });
    }

    updateTransform(state: State, a: Transform.Ref, params: any) {
        const tree = state.build().to(a).update(params);
        return PluginCommands.State.Update.dispatch(this, { state, tree });
    }

    _test_createState(url: string) {
        const b = this.state.data.tree.build();

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
            .apply(StateTransforms.Model.ParseTrajectoryFromMmCif, {}, { ref: 'trajectory' })
            .apply(StateTransforms.Model.CreateModelFromTrajectory, { modelIndex: 0 }, { ref: 'model' })
            .apply(StateTransforms.Model.CreateStructureFromModel, { }, { ref: 'structure' })
            .apply(StateTransforms.Model.CreateStructureAssembly)
            .apply(StateTransforms.Model.CreateStructureSelection, { query, label: 'ALA residues' })
            .apply(StateTransforms.Visuals.CreateStructureRepresentation)
            .getTree();

        this.runTask(this.state.data.update(newTree));
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

    async _test_nextModel() {
        const traj = this.state.data.select('trajectory')[0].obj as SO.Molecule.Trajectory;
        //const modelIndex = (this.state.data.select('model')[0].transform.params as CreateModelFromTrajectory.Params).modelIndex;
        const newTree = this.state.data.build().to('model').update(
            StateTransforms.Model.CreateModelFromTrajectory,
            old => ({ modelIndex: (old.modelIndex + 1) % traj.data.length }))
            .getTree();
        // const newTree = StateTree.updateParams(this.state.data.tree, 'model', { modelIndex: (modelIndex + 1) % traj.data.length });
        await this.runTask(this.state.data.update(newTree));
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