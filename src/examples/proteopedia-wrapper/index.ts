/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createPlugin, DefaultPluginSpec } from 'mol-plugin';
import './index.html'
import { PluginContext } from 'mol-plugin/context';
import { PluginCommands } from 'mol-plugin/command';
import { StateTransforms } from 'mol-plugin/state/transforms';
import { StructureRepresentation3DHelpers } from 'mol-plugin/state/transforms/representation';
import { Color } from 'mol-util/color';
import { PluginStateObject as PSO, PluginStateObject } from 'mol-plugin/state/objects';
import { AnimateModelIndex } from 'mol-plugin/state/animation/built-in';
import { StateBuilder, StateObject } from 'mol-state';
import { EvolutionaryConservation } from './annotation';
import { LoadParams, SupportedFormats, RepresentationStyle, ModelInfo, StateElements } from './helpers';
import { RxEventHelper } from 'mol-util/rx-event-helper';
import { ControlsWrapper } from './ui/controls';
import { PluginState } from 'mol-plugin/state';
import { Scheduler } from 'mol-task';
import { createProteopediaCustomTheme } from './coloring';
require('mol-plugin/skin/light.scss')

class MolStarProteopediaWrapper {
    static VERSION_MAJOR = 2;
    static VERSION_MINOR = 0;

    private _ev = RxEventHelper.create();

    readonly events = {
        modelInfo: this._ev<ModelInfo>()
    };

    plugin: PluginContext;

    init(target: string | HTMLElement, options?: {
        customColorList?: number[],
        customColorDefault?: number
    }) {
        this.plugin = createPlugin(typeof target === 'string' ? document.getElementById(target)! : target, {
            ...DefaultPluginSpec,
            animations: [
                AnimateModelIndex
            ],
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: false
                },
                controls: {
                    right: ControlsWrapper
                }
            }
        });

        const customColoring = createProteopediaCustomTheme((options && options.customColorList) || [], (options && options.customColorDefault) || 0x777777);

        this.plugin.structureRepresentation.themeCtx.colorThemeRegistry.add('proteopedia-custom', customColoring);
        this.plugin.structureRepresentation.themeCtx.colorThemeRegistry.add(EvolutionaryConservation.Descriptor.name, EvolutionaryConservation.colorTheme!);
        this.plugin.lociLabels.addProvider(EvolutionaryConservation.labelProvider);
        this.plugin.customModelProperties.register(EvolutionaryConservation.propertyProvider);
    }

    get state() {
        return this.plugin.state.dataState;
    }

    private download(b: StateBuilder.To<PSO.Root>, url: string) {
        return b.apply(StateTransforms.Data.Download, { url, isBinary: false })
    }

    private model(b: StateBuilder.To<PSO.Data.Binary | PSO.Data.String>, format: SupportedFormats, assemblyId: string) {
        const parsed = format === 'cif'
            ? b.apply(StateTransforms.Data.ParseCif).apply(StateTransforms.Model.TrajectoryFromMmCif)
            : b.apply(StateTransforms.Model.TrajectoryFromPDB);

        return parsed
            .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 }, { ref: StateElements.Model });
    }

    private structure(assemblyId: string) {
        const model = this.state.build().to(StateElements.Model);

        const s = model
            .apply(StateTransforms.Model.CustomModelProperties, { properties: [EvolutionaryConservation.Descriptor.name] }, { ref: StateElements.ModelProps, state: { isGhost: false } })
            .apply(StateTransforms.Model.StructureAssemblyFromModel, { id: assemblyId || 'deposited' }, { ref: StateElements.Assembly });

        s.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-sequence' }, { ref: StateElements.Sequence });
        s.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-het' }, { ref: StateElements.Het });
        s.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' }, { ref: StateElements.Water });

        return s;
    }

    private visual(_style?: RepresentationStyle, partial?: boolean) {
        const structure = this.getObj<PluginStateObject.Molecule.Structure>(StateElements.Assembly);
        if (!structure) return;

        const style = _style || { };

        const update = this.state.build();

        if (!partial || (partial && style.sequence)) {
            const root = update.to(StateElements.Sequence);
            if (style.sequence && style.sequence.hide) {
                root.delete(StateElements.SequenceVisual);
            } else {
                root.applyOrUpdate(StateElements.SequenceVisual, StateTransforms.Representation.StructureRepresentation3D,
                    StructureRepresentation3DHelpers.getDefaultParamsWithTheme(this.plugin,
                        (style.sequence && style.sequence.kind) || 'cartoon',
                        (style.sequence && style.sequence.coloring) || 'unit-index', structure));
            }
        }

        if (!partial || (partial && style.hetGroups)) {
            const root = update.to(StateElements.Het);
            if (style.hetGroups && style.hetGroups.hide) {
                root.delete(StateElements.HetVisual);
            } else {
                if (style.hetGroups && style.hetGroups.hide) {
                    root.delete(StateElements.HetVisual);
                } else {
                    root.applyOrUpdate(StateElements.HetVisual, StateTransforms.Representation.StructureRepresentation3D,
                        StructureRepresentation3DHelpers.getDefaultParamsWithTheme(this.plugin,
                            (style.hetGroups && style.hetGroups.kind) || 'ball-and-stick',
                            (style.hetGroups && style.hetGroups.coloring), structure));
                }
            }
        }

        if (!partial || (partial && style.snfg3d)) {
            const root = update.to(StateElements.Het);
            if (style.hetGroups && style.hetGroups.hide) {
                root.delete(StateElements.HetVisual);
            } else {
                if (style.snfg3d && style.snfg3d.hide) {
                    root.delete(StateElements.Het3DSNFG);
                } else {
                    root.applyOrUpdate(StateElements.Het3DSNFG, StateTransforms.Representation.StructureRepresentation3D,
                        StructureRepresentation3DHelpers.getDefaultParamsWithTheme(this.plugin, 'carbohydrate', void 0, structure));
                }
            }
        }

        if (!partial || (partial && style.water)) {
            const root = update.to(StateElements.Het);
            if (style.water && style.water.hide) {
                root.delete(StateElements.Water);
            } else {
                root.applyOrUpdate(StateElements.Water, StateTransforms.Model.StructureComplexElement, { type: 'water' })
                    .applyOrUpdate(StateElements.WaterVisual, StateTransforms.Representation.StructureRepresentation3D,
                        StructureRepresentation3DHelpers.getDefaultParamsWithTheme(this.plugin,
                            (style.water && style.water.kind) || 'ball-and-stick',
                            (style.water && style.water.coloring), structure, { alpha: 0.51 }));
            }
        }

        return update;
    }

    private getObj<T extends StateObject>(ref: string): T['data'] {
        const state = this.state;
        const cell = state.select(ref)[0];
        if (!cell || !cell.obj) return void 0;
        return (cell.obj as T).data;
    }

    private async doInfo(checkPreferredAssembly: boolean) {
        const model = this.getObj<PluginStateObject.Molecule.Model>('model');
        if (!model) return;

        const info = await ModelInfo.get(this.plugin, model, checkPreferredAssembly)
        this.events.modelInfo.next(info);
        return info;
    }

    private applyState(tree: StateBuilder) {
        return PluginCommands.State.Update.dispatch(this.plugin, { state: this.plugin.state.dataState, tree });
    }

    private loadedParams: LoadParams = { url: '', format: 'cif', assemblyId: '' };
    async load({ url, format = 'cif', assemblyId = '', representationStyle }: LoadParams) {
        let loadType: 'full' | 'update' = 'full';

        const state = this.plugin.state.dataState;

        if (this.loadedParams.url !== url || this.loadedParams.format !== format) {
            loadType = 'full';
        } else if (this.loadedParams.url === url) {
            if (state.select(StateElements.Assembly).length > 0) loadType = 'update';
        }

        if (loadType === 'full') {
            await PluginCommands.State.RemoveObject.dispatch(this.plugin, { state, ref: state.tree.root.ref });
            const modelTree = this.model(this.download(state.build().toRoot(), url), format, assemblyId);
            await this.applyState(modelTree);
            const info = await this.doInfo(true);
            const structureTree = this.structure((assemblyId === 'preferred' && info && info.preferredAssemblyId) || assemblyId);
            await this.applyState(structureTree);
        } else {
            const tree = state.build();
            tree.to(StateElements.Assembly).update(StateTransforms.Model.StructureAssemblyFromModel, p => ({ ...p, id: assemblyId || 'deposited' }));
            await this.applyState(tree);
        }

        await this.updateStyle(representationStyle);

        this.loadedParams = { url, format, assemblyId };
        Scheduler.setImmediate(() => PluginCommands.Camera.Reset.dispatch(this.plugin, { }));
    }

    async updateStyle(style?: RepresentationStyle, partial?: boolean) {
        const tree = this.visual(style, partial);
        if (!tree) return;
        await PluginCommands.State.Update.dispatch(this.plugin, { state: this.plugin.state.dataState, tree });
    }

    setBackground(color: number) {
        PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: { backgroundColor: Color(color) } });
    }

    toggleSpin() {
        const trackball = this.plugin.canvas3d.props.trackball;
        const spinning = trackball.spin;
        PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: { trackball: { ...trackball, spin: !trackball.spin } } });
        if (!spinning) PluginCommands.Camera.Reset.dispatch(this.plugin, { });
    }

    animate = {
        modelIndex: {
            maxFPS: 8,
            onceForward: () => { this.plugin.state.animation.play(AnimateModelIndex, { maxFPS: Math.max(0.5, this.animate.modelIndex.maxFPS | 0), mode: { name: 'once', params: { direction: 'forward' } } }) },
            onceBackward: () => { this.plugin.state.animation.play(AnimateModelIndex, { maxFPS: Math.max(0.5, this.animate.modelIndex.maxFPS | 0), mode: { name: 'once', params: { direction: 'backward' } } }) },
            palindrome: () => { this.plugin.state.animation.play(AnimateModelIndex, { maxFPS: Math.max(0.5, this.animate.modelIndex.maxFPS | 0), mode: { name: 'palindrome', params: {} } }) },
            loop: () => { this.plugin.state.animation.play(AnimateModelIndex, { maxFPS: Math.max(0.5, this.animate.modelIndex.maxFPS | 0), mode: { name: 'loop', params: {} } }) },
            stop: () => this.plugin.state.animation.stop()
        }
    }

    coloring = {
        evolutionaryConservation: async () => {
            await this.updateStyle({ sequence: { kind: 'spacefill' } }, true);

            const state = this.state;

            // const visuals = state.selectQ(q => q.ofType(PluginStateObject.Molecule.Structure.Representation3D).filter(c => c.transform.transformer === StateTransforms.Representation.StructureRepresentation3D));
            const tree = state.build();
            const colorTheme = { name: EvolutionaryConservation.Descriptor.name, params: this.plugin.structureRepresentation.themeCtx.colorThemeRegistry.get(EvolutionaryConservation.Descriptor.name).defaultValues };

            tree.to(StateElements.SequenceVisual).update(StateTransforms.Representation.StructureRepresentation3D, old => ({ ...old, colorTheme }));
            // for (const v of visuals) {
            // }

            await PluginCommands.State.Update.dispatch(this.plugin, { state, tree });
        }
    }

    snapshot = {
        get: () => {
            return this.plugin.state.getSnapshot();
        },
        set: (snapshot: PluginState.Snapshot) => {
            return this.plugin.state.setSnapshot(snapshot);
        },
        download: async (url: string) => {
            try {
                const data = await this.plugin.runTask(this.plugin.fetch({ url }));
                const snapshot = JSON.parse(data);
                await this.plugin.state.setSnapshot(snapshot);
            } catch (e) {
                console.log(e);
            }
        }

    }
}

(window as any).MolStarProteopediaWrapper = MolStarProteopediaWrapper;