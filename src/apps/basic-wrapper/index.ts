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
import { StateTree } from 'mol-state';
import { Color } from 'mol-util/color';
require('mol-plugin/skin/light.scss')

class BasicWrapper {
    plugin: PluginContext;
    stateTemplate: StateTree;

    init(target: string | HTMLElement) {
        this.plugin = createPlugin(typeof target === 'string' ? document.getElementById(target)! : target, {
            ...DefaultPluginSpec,
            initialLayout: {
                isExpanded: false,
                showControls: false
            }
        });

        const state = this.plugin.state.dataState.build();
        const visualRoot = state.toRoot()
            .apply(StateTransforms.Data.Download, { url: '', isBinary: false }, { ref: 'url' })
            .apply(StateTransforms.Data.ParseCif)
            .apply(StateTransforms.Model.TrajectoryFromMmCif)
            .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 })
            .apply(StateTransforms.Model.StructureAssemblyFromModel, { id: '' }, { ref: 'asm' })

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

        this.stateTemplate = state.getTree();
    }

    async loadCif(url: string, assemblyId?: string) {
        const state = this.stateTemplate.build();

        state.to('url').update(StateTransforms.Data.Download, p => ({ ...p, url }));
        state.to('asm').update(StateTransforms.Model.StructureAssemblyFromModel, p => ({ ...p, id: assemblyId }));

        await PluginCommands.State.Update.dispatch(this.plugin, { state: this.plugin.state.dataState, tree: state });

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
}

(window as any).BasicMolStarWrapper = new BasicWrapper();