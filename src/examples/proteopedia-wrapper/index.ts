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
import { LoadParams, SupportedFormats, RepresentationStyle, ModelInfo } from './helpers';
import { RxEventHelper } from 'mol-util/rx-event-helper';
require('mol-plugin/skin/light.scss')

class MolStarProteopediaWrapper {
    static VERSION_MAJOR = 1;

    private _ev = RxEventHelper.create();

    readonly events = {
        modelInfo: this._ev<ModelInfo>()
    };

    plugin: PluginContext;

    init(target: string | HTMLElement) {
        this.plugin = createPlugin(typeof target === 'string' ? document.getElementById(target)! : target, {
            ...DefaultPluginSpec,
            initialLayout: {
                isExpanded: false,
                showControls: false
            }
        });

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
            .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 }, { ref: 'model' });
    }

    private structure(assemblyId: string) {
        const model = this.state.build().to('model');

        return model
            .apply(StateTransforms.Model.CustomModelProperties, { properties: [EvolutionaryConservation.Descriptor.name] }, { ref: 'props', props: { isGhost: false } })
            .apply(StateTransforms.Model.StructureAssemblyFromModel, { id: assemblyId || 'deposited' }, { ref: 'asm' });
    }

    private visual(ref: string, style?: RepresentationStyle) {
        const structure = this.getObj<PluginStateObject.Molecule.Structure>(ref);
        if (!structure) return;

        const root = this.state.build().to(ref);

        root.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-sequence' }, { ref: 'sequence' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsWithTheme(this.plugin,
                    (style && style.sequence && style.sequence.kind) || 'cartoon',
                    (style && style.sequence && style.sequence.coloring) || 'unit-index', structure),
                    { ref: 'sequence-visual' });
        root.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-het' }, { ref: 'het' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsWithTheme(this.plugin,
                    (style && style.hetGroups && style.hetGroups.kind) || 'ball-and-stick',
                    (style && style.hetGroups && style.hetGroups.coloring), structure),
                    { ref: 'het-visual' });
        root.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' }, { ref: 'water' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsWithTheme(this.plugin,
                    (style && style.water && style.water.kind) || 'ball-and-stick',
                    (style && style.water && style.water.coloring), structure, { alpha: 0.51 }),
                    { ref: 'water-visual' });

        return root;
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
            if (state.select('asm').length > 0) loadType = 'update';
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
            tree.to('asm').update(StateTransforms.Model.StructureAssemblyFromModel, p => ({ ...p, id: assemblyId || 'deposited' }));
            await this.applyState(tree);
        }

        await this.updateStyle(representationStyle);

        this.loadedParams = { url, format, assemblyId };
        PluginCommands.Camera.Reset.dispatch(this.plugin, { });
    }

    async updateStyle(style?: RepresentationStyle) {
        const tree = this.visual('asm', style);
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
            await this.updateStyle({ sequence: { kind: 'molecular-surface' } });

            const state = this.state;

            //const visuals = state.selectQ(q => q.ofType(PluginStateObject.Molecule.Representation3D).filter(c => c.transform.transformer === StateTransforms.Representation.StructureRepresentation3D));
            const tree = state.build();
            const colorTheme = { name: EvolutionaryConservation.Descriptor.name, params: this.plugin.structureRepresentation.themeCtx.colorThemeRegistry.get(EvolutionaryConservation.Descriptor.name).defaultValues };

            tree.to('sequence-visual').update(StateTransforms.Representation.StructureRepresentation3D, old => ({ ...old, colorTheme }));
            //for (const v of visuals) {
            //}

            await PluginCommands.State.Update.dispatch(this.plugin, { state, tree });
        }
    }
}

(window as any).MolStarProteopediaWrapper = MolStarProteopediaWrapper;