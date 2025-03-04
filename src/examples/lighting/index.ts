/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Canvas3DProps } from '../../mol-canvas3d/canvas3d';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { createPluginUI } from '../../mol-plugin-ui';
import { PluginUIContext } from '../../mol-plugin-ui/context';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginCommands } from '../../mol-plugin/commands';
import { Asset } from '../../mol-util/assets';
import { Color } from '../../mol-util/color';
import './index.html';
import '../../mol-plugin-ui/skin/light.scss';

type LoadParams = { url: string, format?: BuiltInTrajectoryFormat, isBinary?: boolean, assemblyId?: string }

type _Preset = Pick<Canvas3DProps, 'postprocessing' | 'renderer'>
type Preset = { [K in keyof _Preset]: Partial<_Preset[K]> }

const Canvas3DPresets = {
    illustrative: {
        canvas3d: <Preset>{
            postprocessing: {
                occlusion: {
                    name: 'on',
                    params: {
                        samples: 32,
                        multiScale: { name: 'off', params: {} },
                        radius: 5,
                        bias: 0.8,
                        blurKernelSize: 15,
                        blurDepthBias: 0.5,
                        resolutionScale: 1,
                        color: Color(0x000000),
                        transparentThreshold: 0.4,
                    }
                },
                outline: {
                    name: 'on',
                    params: {
                        scale: 1,
                        threshold: 0.33,
                        color: Color(0x000000),
                        includeTransparent: true,
                    }
                },
                shadow: {
                    name: 'off',
                    params: {}
                },
            },
            renderer: {
                ambientIntensity: 1.0,
                light: []
            }
        }
    },
    occlusion: {
        canvas3d: <Preset>{
            postprocessing: {
                occlusion: {
                    name: 'on',
                    params: {
                        samples: 32,
                        multiScale: { name: 'off', params: {} },
                        radius: 5,
                        bias: 0.8,
                        blurKernelSize: 15,
                        resolutionScale: 1,
                    }
                },
                outline: {
                    name: 'off',
                    params: {}
                },
                shadow: {
                    name: 'off',
                    params: {}
                },
            },
            renderer: {
                ambientIntensity: 0.4,
                light: [{ inclination: 180, azimuth: 0, color: Color.fromNormalizedRgb(1.0, 1.0, 1.0),
                    intensity: 0.6 }]
            }
        }
    },
    standard: {
        canvas3d: <Preset>{
            postprocessing: {
                occlusion: { name: 'off', params: {} },
                outline: { name: 'off', params: {} },
                shadow: { name: 'off', params: {} },
            },
            renderer: {
                ambientIntensity: 0.4,
                light: [{ inclination: 180, azimuth: 0, color: Color.fromNormalizedRgb(1.0, 1.0, 1.0),
                    intensity: 0.6 }]
            }
        }
    }
};


type Canvas3DPreset = keyof typeof Canvas3DPresets

class LightingDemo {
    plugin: PluginUIContext;

    private radius = 5;
    private bias = 1.1;
    private preset: Canvas3DPreset = 'illustrative';

    async init(target: string | HTMLElement) {
        this.plugin = await createPluginUI({
            target: typeof target === 'string' ? document.getElementById(target)! : target,
            render: renderReact18,
            spec: {
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
            }
        });

        this.setPreset('illustrative');
    }

    setPreset(preset: Canvas3DPreset) {
        const props = Canvas3DPresets[preset];
        if (props.canvas3d.postprocessing.occlusion?.name === 'on') {
            props.canvas3d.postprocessing.occlusion.params.radius = this.radius;
            props.canvas3d.postprocessing.occlusion.params.bias = this.bias;
        }
        PluginCommands.Canvas3D.SetSettings(this.plugin, {
            settings: {
                ...props,
                renderer: {
                    ...this.plugin.canvas3d!.props.renderer,
                    ...props.canvas3d.renderer
                },
                postprocessing: {
                    ...this.plugin.canvas3d!.props.postprocessing,
                    ...props.canvas3d.postprocessing
                },
            }
        });
    }

    async load({ url, format = 'mmcif', isBinary = true, assemblyId = '' }: LoadParams, radius: number, bias: number) {
        await this.plugin.clear();

        const data = await this.plugin.builders.data.download({ url: Asset.Url(url), isBinary }, { state: { isGhost: true } });
        const trajectory = await this.plugin.builders.structure.parseTrajectory(data, format);
        const model = await this.plugin.builders.structure.createModel(trajectory);
        const structure = await this.plugin.builders.structure.createStructure(model, assemblyId ? { name: 'assembly', params: { id: assemblyId } } : { name: 'model', params: {} });

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