/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Canvas3DProps } from '../../mol-canvas3d/canvas3d';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { createPluginUI } from '../../mol-plugin-ui/react18';
import { PluginUIContext } from '../../mol-plugin-ui/context';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginCommands } from '../../mol-plugin/commands';
import { Asset } from '../../mol-util/assets';
import { Color } from '../../mol-util/color';
import './index.html';
import { Structure, StructureElement, StructureProperties } from '../../mol-model/structure';
import { StructureMeasurementManager } from '../../mol-plugin-state/manager/structure/measurement';
import { Loci } from '../../mol-model/structure/structure/element/loci';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { SequenceWrapper } from '../../mol-plugin-ui/sequence/wrapper';
import { getChainOptions, getModelEntityOptions, getOperatorOptions, getSequenceWrapper } from '../../mol-plugin-ui/sequence';
require('mol-plugin-ui/skin/light.scss');

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
                        resolutionScale: 1,
                        color: Color(0x000000),
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
        this.plugin = await createPluginUI(typeof target === 'string' ? document.getElementById(target)! : target, {
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
        let structure: any = await this.plugin.builders.structure.createStructure(model, assemblyId ? { name: 'assembly', params: { id: assemblyId } } : { name: 'model', params: {} });

        const polymer = await this.plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer');
        if (polymer) await this.plugin.builders.structure.representation.addRepresentation(polymer, { type: 'cartoon', color: 'illustrative' });

        const atomLoci: any = [];
        Structure.eachAtomicHierarchyElement(this.plugin.managers.structure.hierarchy.selection.structures[0].cell.obj?.data as Structure, {
            atom: (a) => {
                const newloci = Structure.toStructureElementLoci(a.structure);
                console.log('newLoci', newloci);
                atomLoci.push(newloci);

                console.log(StructureProperties.atom.x(a), StructureProperties.atom.y(a), StructureProperties.atom.z(a));
            },
        });

        setTimeout(() => {
            const state = this.plugin.state.data;
            const cell = state.select(structure.ref)[0];
            if (!structure.ref || !cell || !cell.obj) return;
            structure = (cell.obj as PluginStateObject.Molecule.Structure).data;
            console.log('structure', structure)
            const wrappers: { wrapper: (string | SequenceWrapper.Any), label: string }[] = [];
            // Get wrapper of structure

            for (const [modelEntityId, eLabel] of getModelEntityOptions(structure, true)) {
                for (const [chainGroupId, cLabel] of getChainOptions(structure, modelEntityId)) {
                    for (const [operatorKey] of getOperatorOptions(structure, modelEntityId, chainGroupId)) {
                        wrappers.push({
                            wrapper: getSequenceWrapper({
                                structure,
                                modelEntityId,
                                chainGroupId,
                                operatorKey
                            }, this.plugin.managers.structure.selection),
                            label: `${cLabel} | ${eLabel}`
                        });
                        if (wrappers.length > 30) return [];
                    }
                }
            }

            const wrapper = wrappers[0];
            const sequenceWrapper = wrapper.wrapper as SequenceWrapper.Any;
            const measurementManager = new StructureMeasurementManager(this.plugin);
            console.log('Wrapper', sequenceWrapper.length);
            measurementManager.addDistance(sequenceWrapper.getLoci(80), sequenceWrapper.getLoci(100));
            // const loci = sequenceWrapper.getLoci(i);
        }, 3000);

        this.radius = radius;
        this.bias = bias;
        this.setPreset(this.preset);
    }

}

function toEntries<T>(a: T[]) {
    return a.map((value, index) => [index, value] as const);
}
(window as any).LightingDemo = new LightingDemo();