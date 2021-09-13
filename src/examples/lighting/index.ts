/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Canvas3DProps } from '../../mol-canvas3d/canvas3d';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { createPlugin } from '../../mol-plugin-ui';
import { PluginUIContext } from '../../mol-plugin-ui/context';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginCommands } from '../../mol-plugin/commands';
import { Asset } from '../../mol-util/assets';
import './index.html';
require('mol-plugin-ui/skin/light.scss');

type LoadParams = { url: string, format?: BuiltInTrajectoryFormat, isBinary?: boolean, assemblyId?: string }

type _Preset = Pick<Canvas3DProps, 'multiSample' | 'postprocessing' | 'renderer'>
type Preset = { [K in keyof _Preset]: Partial<_Preset[K]> }

const Canvas3DPresets = {
    illustrative: <Preset> {
        multiSample: {
            mode: 'temporal' as Canvas3DProps['multiSample']['mode']
        },
        postprocessing: {
            occlusion: { name: 'on', params: { samples: 32, radius: 6, bias: 1.4, blurKernelSize: 15 } },
            outline: { name: 'on', params: { scale: 1, threshold: 0.1 } }
        },
        renderer: {
            style: { name: 'flat', params: {} }
        }
    },
    occlusion: <Preset> {
        multiSample: {
            mode: 'temporal' as Canvas3DProps['multiSample']['mode']
        },
        postprocessing: {
            occlusion: { name: 'on', params: { samples: 32, radius: 6, bias: 1.4, blurKernelSize: 15 } },
            outline: { name: 'off', params: { } }
        },
        renderer: {
            style: { name: 'matte', params: {} }
        }
    },
    standard: <Preset> {
        multiSample: {
            mode: 'off' as Canvas3DProps['multiSample']['mode']
        },
        postprocessing: {
            occlusion: { name: 'off', params: { } },
            outline: { name: 'off', params: { } }
        },
        renderer: {
            style: { name: 'matte', params: {} }
        }
    }
};


type Canvas3DPreset = keyof typeof Canvas3DPresets

class LightingDemo {
    plugin: PluginUIContext;

    private radius = 5;
    private bias = 1.1;
    private preset: Canvas3DPreset = 'illustrative';

    init(target: string | HTMLElement) {
        this.plugin = createPlugin(typeof target === 'string' ? document.getElementById(target)! : target, {
            ...DefaultPluginUISpec(),
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: false
                },
            },
            components: {
                controls: { left: 'none', right: 'none', top: 'none', bottom: 'none' }
            }
        });

        this.setPreset('illustrative');
    }

    setPreset(preset: Canvas3DPreset) {
        const props = Canvas3DPresets[preset];
        if (props.postprocessing.occlusion?.name === 'on') {
            props.postprocessing.occlusion.params.radius = this.radius;
            props.postprocessing.occlusion.params.bias = this.bias;
        }
        PluginCommands.Canvas3D.SetSettings(this.plugin, { settings: {
            ...props,
            multiSample: {
                ...this.plugin.canvas3d!.props.multiSample,
                ...props.multiSample
            },
            renderer: {
                ...this.plugin.canvas3d!.props.renderer,
                ...props.renderer
            },
            postprocessing: {
                ...this.plugin.canvas3d!.props.postprocessing,
                ...props.postprocessing
            },
        } });
    }

    async load({ url, format = 'mmcif', isBinary = true, assemblyId = '' }: LoadParams, radius: number, bias: number) {
        await this.plugin.clear();

        const data = await this.plugin.builders.data.download({ url: Asset.Url(url), isBinary }, { state: { isGhost: true } });
        const trajectory = await this.plugin.builders.structure.parseTrajectory(data, format);
        const model = await this.plugin.builders.structure.createModel(trajectory);
        const structure = await this.plugin.builders.structure.createStructure(model, assemblyId ? { name: 'assembly', params: { id: assemblyId } } : { name: 'model', params: { } });

        const polymer = await this.plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer');
        if (polymer) await this.plugin.builders.structure.representation.addRepresentation(polymer, { type: 'spacefill', color: 'illustrative' });

        const ligand = await this.plugin.builders.structure.tryCreateComponentStatic(structure, 'ligand');
        if (ligand) await this.plugin.builders.structure.representation.addRepresentation(ligand, { type: 'ball-and-stick', color: 'element-symbol', colorParams: { carbonColor: { name: 'element-symbol', params: {} } } });

        this.radius = radius;
        this.bias = bias;
        this.setPreset(this.preset);
    }
}

(window as any).LightingDemo = new LightingDemo();