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
import { StateTreeBuilder } from 'mol-state/tree/builder';
import { PluginStateObject as PSO } from 'mol-plugin/state/objects';
import { AnimateModelIndex } from 'mol-plugin/state/animation/built-in';
require('mol-plugin/skin/light.scss')

type SupportedFormats = 'cif' | 'pdb'
type LoadParams = { url: string, format?: SupportedFormats, assemblyId?: string }

class BasicWrapper {
    plugin: PluginContext;

    init(target: string | HTMLElement) {
        this.plugin = createPlugin(typeof target === 'string' ? document.getElementById(target)! : target, {
            ...DefaultPluginSpec,
            initialLayout: {
                isExpanded: false,
                showControls: false
            }
        });
    }

    private download(b: StateTreeBuilder.To<PSO.Root>, url: string) {
        return b.apply(StateTransforms.Data.Download, { url, isBinary: false })
    }

    private parse(b: StateTreeBuilder.To<PSO.Data.Binary | PSO.Data.String>, format: SupportedFormats, assemblyId: string) {
        const parsed = format === 'cif'
            ? b.apply(StateTransforms.Data.ParseCif).apply(StateTransforms.Model.TrajectoryFromMmCif)
            : b.apply(StateTransforms.Model.TrajectoryFromPDB);

        return parsed
            .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 })
            .apply(StateTransforms.Model.StructureAssemblyFromModel, { id: assemblyId || 'deposited' }, { ref: 'asm' });
    }

    private visual(visualRoot: StateTreeBuilder.To<PSO.Molecule.Structure>) {
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-sequence' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(this.plugin, 'cartoon'));
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-het' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(this.plugin, 'ball-and-stick'));
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(this.plugin, 'ball-and-stick', { alpha: 0.51 }));
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'spheres' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(this.plugin, 'spacefill'));
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

        let tree: StateTreeBuilder.Root;
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
            forward: (maxFPS = 8) => { this.plugin.state.animation.play(AnimateModelIndex, { maxFPS: Math.max(0.5, maxFPS | 0), direction: 'forward' }) },
            backward: (maxFPS = 8) => { this.plugin.state.animation.play(AnimateModelIndex, { maxFPS: Math.max(0.5, maxFPS | 0), direction: 'backward' }) },
            stop: () => this.plugin.state.animation.stop()
        }
    }
}

(window as any).BasicMolStarWrapper = new BasicWrapper();