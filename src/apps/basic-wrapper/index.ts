/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createPlugin, DefaultPluginSpec } from '../../mol-plugin';
import './index.html'
import { PluginContext } from '../../mol-plugin/context';
import { PluginCommands } from '../../mol-plugin/command';
import { StateTransforms } from '../../mol-plugin/state/transforms';
import { StructureRepresentation3DHelpers } from '../../mol-plugin/state/transforms/representation';
import { Color } from '../../mol-util/color';
import { PluginStateObject as PSO, PluginStateObject } from '../../mol-plugin/state/objects';
import { AnimateModelIndex } from '../../mol-plugin/state/animation/built-in';
import { StateBuilder, StateTransform } from '../../mol-state';
import { StripedResidues } from './coloring';
import { StaticSuperpositionTestData, buildStaticSuperposition, dynamicSuperpositionTest } from './superposition';
import { PDBeStructureQualityReport } from '../../mol-plugin/behavior/dynamic/custom-props';
import { CustomToastMessage } from './controls';
import { EmptyLoci } from '../../mol-model/loci';
import { StructureSelection } from '../../mol-model/structure';
import { Script } from '../../mol-script/script';
require('mol-plugin/skin/light.scss')

type SupportedFormats = 'cif' | 'pdb'
type LoadParams = { url: string, format?: SupportedFormats, assemblyId?: string }

class BasicWrapper {
    plugin: PluginContext;

    init(target: string | HTMLElement) {
        this.plugin = createPlugin(typeof target === 'string' ? document.getElementById(target)! : target, {
            ...DefaultPluginSpec,
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: false
                },
                controls: {
                    // left: 'none',
                    // right: BasicWrapperControls
                }
            },
            components: {
                remoteState: 'none'
            }
        });

        this.plugin.structureRepresentation.themeCtx.colorThemeRegistry.add(StripedResidues.Descriptor.name, StripedResidues.colorTheme!);
        this.plugin.lociLabels.addProvider(StripedResidues.labelProvider);
        this.plugin.customModelProperties.register(StripedResidues.propertyProvider);
    }

    private download(b: StateBuilder.To<PSO.Root>, url: string) {
        return b.apply(StateTransforms.Data.Download, { url, isBinary: false })
    }

    private parse(b: StateBuilder.To<PSO.Data.Binary | PSO.Data.String>, format: SupportedFormats, assemblyId: string) {
        const parsed = format === 'cif'
            ? b.apply(StateTransforms.Data.ParseCif).apply(StateTransforms.Model.TrajectoryFromMmCif)
            : b.apply(StateTransforms.Model.TrajectoryFromPDB);

        return parsed
            .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 })
            .apply(StateTransforms.Model.CustomModelProperties, { properties: [StripedResidues.Descriptor.name] }, { ref: 'props', state: { isGhost: false } })
            .apply(StateTransforms.Model.StructureAssemblyFromModel, { id: assemblyId || 'deposited' }, { ref: 'asm' });
    }

    private visual(visualRoot: StateBuilder.To<PSO.Molecule.Structure>) {
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-sequence' }, { ref: 'seq' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(this.plugin, 'cartoon'), { ref: 'seq-visual' });
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-het' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(this.plugin, 'ball-and-stick'), { ref: 'het-visual' });
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(this.plugin, 'ball-and-stick', { alpha: 0.51 }), { ref: 'water-visual' });
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'spheres' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(this.plugin, 'spacefill'), { ref: 'ihm-visual' });
        return visualRoot;
    }

    private loadedParams: LoadParams = { url: '', format: 'cif', assemblyId: '' };
    async load({ url, format = 'cif', assemblyId = '' }: LoadParams) {
        let loadType: 'full' | 'update' = 'full';

        const state = this.plugin.state.dataState;

        if (this.loadedParams.url !== url || this.loadedParams.format !== format) {
            loadType = 'full';
        } else if (this.loadedParams.url === url) {
            if (state.select('asm').length > 0) loadType = 'update';
        }

        let tree: StateBuilder.Root;
        if (loadType === 'full') {
            await PluginCommands.State.RemoveObject.dispatch(this.plugin, { state, ref: state.tree.root.ref });
            tree = state.build();
            this.visual(this.parse(this.download(tree.toRoot(), url), format, assemblyId));
        } else {
            tree = state.build();
            tree.to('asm').update(StateTransforms.Model.StructureAssemblyFromModel, p => ({ ...p, id: assemblyId || 'deposited' }));
        }

        await PluginCommands.State.Update.dispatch(this.plugin, { state: this.plugin.state.dataState, tree });
        this.loadedParams = { url, format, assemblyId };
        PluginCommands.Camera.Reset.dispatch(this.plugin, { });
    }

    setBackground(color: number) {
        const renderer = this.plugin.canvas3d!.props.renderer;
        PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: { renderer: { ...renderer,  backgroundColor: Color(color) } } });
    }

    toggleSpin() {
        if (!this.plugin.canvas3d) return;

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
        applyStripes: async () => {
            const state = this.plugin.state.dataState;

            const visuals = state.selectQ(q => q.ofTransformer(StateTransforms.Representation.StructureRepresentation3D));
            const tree = state.build();
            const colorTheme = { name: StripedResidues.Descriptor.name, params: this.plugin.structureRepresentation.themeCtx.colorThemeRegistry.get(StripedResidues.Descriptor.name).defaultValues };

            for (const v of visuals) {
                tree.to(v).update(old => ({ ...old, colorTheme }));
            }

            await PluginCommands.State.Update.dispatch(this.plugin, { state, tree });
        }
    }

    interactivity = {
        highlightOn: () => {
            const seq_id = 7;
            const data = (this.plugin.state.dataState.select('asm')[0].obj as PluginStateObject.Molecule.Structure).data;
            const sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
                'residue-test': Q.core.rel.eq([Q.struct.atomProperty.macromolecular.label_seq_id(), seq_id]),
                'group-by': Q.struct.atomProperty.macromolecular.residueKey()
            }), data);
            const loci = StructureSelection.toLociWithSourceUnits(sel);
            this.plugin.interactivity.lociHighlights.highlightOnly({ loci });
        },
        clearHighlight: () => {
            this.plugin.interactivity.lociHighlights.highlightOnly({ loci: EmptyLoci });
        }
    }

    tests = {
        staticSuperposition: async () => {
            const state = this.plugin.state.dataState;
            const tree = buildStaticSuperposition(this.plugin, StaticSuperpositionTestData);
            await PluginCommands.State.RemoveObject.dispatch(this.plugin, { state, ref: StateTransform.RootRef });
            await PluginCommands.State.Update.dispatch(this.plugin, { state, tree });
        },
        dynamicSuperposition: async () => {
            await PluginCommands.State.RemoveObject.dispatch(this.plugin, { state: this.plugin.state.dataState, ref: StateTransform.RootRef });
            await dynamicSuperpositionTest(this.plugin, ['1tqn', '2hhb', '4hhb'], 'HEM');
        },
        toggleValidationTooltip: async () => {
            const state = this.plugin.state.behaviorState;
            const tree = state.build().to(PDBeStructureQualityReport.id).update(PDBeStructureQualityReport, p => ({ ...p, showTooltip: !p.showTooltip }));
            await PluginCommands.State.Update.dispatch(this.plugin, { state, tree });
        },
        showToasts: () => {
            PluginCommands.Toast.Show.dispatch(this.plugin, {
                title: 'Toast 1',
                message: 'This is an example text, timeout 3s',
                key: 'toast-1',
                timeoutMs: 3000
            });
            PluginCommands.Toast.Show.dispatch(this.plugin, {
                title: 'Toast 2',
                message: CustomToastMessage,
                key: 'toast-2'
            });
        },
        hideToasts: () => {
            PluginCommands.Toast.Hide.dispatch(this.plugin, { key: 'toast-1' });
            PluginCommands.Toast.Hide.dispatch(this.plugin, { key: 'toast-2' });
        }
    }
}

(window as any).BasicMolStarWrapper = new BasicWrapper();