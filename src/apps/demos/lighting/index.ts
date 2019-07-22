/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createPlugin, DefaultPluginSpec } from '../../../mol-plugin';
import './index.html'
import { PluginContext } from '../../../mol-plugin/context';
import { PluginCommands } from '../../../mol-plugin/command';
import { StateTransforms } from '../../../mol-plugin/state/transforms';
import { StructureRepresentation3DHelpers } from '../../../mol-plugin/state/transforms/representation';
import { PluginStateObject as PSO } from '../../../mol-plugin/state/objects';
import { StateBuilder } from '../../../mol-state';
import { Canvas3DProps } from '../../../mol-canvas3d/canvas3d';
require('mol-plugin/skin/light.scss')

type SupportedFormats = 'cif' | 'pdb'
type LoadParams = { url: string, format?: SupportedFormats, assemblyId?: string }

const Canvas3DPresets = {
    illustrative: {
        multiSample: {
            mode: 'temporal' as Canvas3DProps['multiSample']['mode']
        },
        postprocessing: {
            occlusionEnable: true,
            occlusionBias: 0.8,
            occlusionKernelSize: 6,
            outlineEnable: true,
        },
        renderer: {
            ambientIntensity: 1,
            lightIntensity: 0,
        }
    },
    occlusion: {
        multiSample: {
            mode: 'temporal' as Canvas3DProps['multiSample']['mode']
        },
        postprocessing: {
            occlusionEnable: true,
            occlusionBias: 0.8,
            occlusionKernelSize: 6,
            outlineEnable: false,
        },
        renderer: {
            ambientIntensity: 0.4,
            lightIntensity: 0.6,
        }
    },
    standard: {
        multiSample: {
            mode: 'off' as Canvas3DProps['multiSample']['mode']
        },
        postprocessing: {
            occlusionEnable: false,
            outlineEnable: false,
        },
        renderer: {
            ambientIntensity: 0.4,
            lightIntensity: 0.6,
        }
    }
}

type Canvas3DPreset = keyof typeof Canvas3DPresets

function getPreset(preset: Canvas3DPreset) {
    switch (preset) {
        case 'illustrative': return Canvas3DPresets['illustrative']
        case 'standard': return Canvas3DPresets['standard']
        case 'occlusion': return Canvas3DPresets['occlusion']
    }
}

class LightingDemo {
    plugin: PluginContext;

    init(target: string | HTMLElement) {
        this.plugin = createPlugin(typeof target === 'string' ? document.getElementById(target)! : target, {
            ...DefaultPluginSpec,
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: false
                },
                controls: { left: 'none', right: 'none', top: 'none', bottom: 'none' }
            }
        });

        this.setPreset('illustrative');
    }

    setPreset(preset: Canvas3DPreset) {
        const props = getPreset(preset)
        PluginCommands.Canvas3D.SetSettings.dispatch(this.plugin, { settings: {
            ...props,
            multiSample: {
                ...this.plugin.canvas3d.props.multiSample,
                ...props.multiSample
            },
            renderer: {
                ...this.plugin.canvas3d.props.renderer,
                ...props.renderer
            },
            postprocessing: {
                ...this.plugin.canvas3d.props.postprocessing,
                ...props.postprocessing
            },
        }});
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
            .apply(StateTransforms.Model.StructureAssemblyFromModel, { id: assemblyId || 'deposited' }, { ref: 'asm' });
    }

    private visual(visualRoot: StateBuilder.To<PSO.Molecule.Structure>) {
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-sequence' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(this.plugin, 'spacefill', {}, 'illustrative'), { ref: 'seq-visual' });
        visualRoot.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-het' })
            .apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsStatic(this.plugin, 'ball-and-stick'), { ref: 'het-visual' });
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
}

(window as any).LightingDemo = new LightingDemo();