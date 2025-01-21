/**
 * Copyright (c) 2022-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PostprocessingParams } from '../../mol-canvas3d/passes/postprocessing';
import { PresetStructureRepresentations } from '../../mol-plugin-state/builder/structure/representation-preset';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { CollapsableControls, PurePluginUIComponent } from '../base';
import { Button } from '../controls/common';
import { MagicWandSvg } from '../controls/icons';


export class StructureQuickStylesControls extends CollapsableControls {
    defaultState() {
        return {
            isCollapsed: false,
            header: 'Quick Styles',
            brand: { accent: 'gray' as const, svg: MagicWandSvg }
        };
    }

    renderControls() {
        return <>
            <QuickStyles />
        </>;
    }
}


type PresetName = 'default' | 'cartoon' | 'spacefill' | 'surface';
type StyleName = 'default' | 'illustrative';

interface QuickStylesState {
    busy: boolean,
    style: StyleName,
}


export class QuickStyles extends PurePluginUIComponent<{}, QuickStylesState> {
    state: QuickStylesState = { busy: false, style: 'default' };

    async applyRepresentation(preset: PresetName) {
        this.setState({ busy: true });
        await applyRepresentationPreset(this.plugin, preset);
        await applyStyle(this.plugin, this.state.style); // reapplying current style is desired because some presets come with weird params (namely spacefill comes with ignoreLight:true)
        this.setState({ busy: false });
    }

    async applyStyle(style: StyleName) {
        this.setState({ busy: true });
        await applyStyle(this.plugin, style);
        this.setState({ busy: false, style });
    }

    render() {
        return <>
            <NoncollapsableGroup title='Apply Representation'>
                <div className='msp-flex-row'>
                    <Button title='Applies default representation preset (depends on structure size)'
                        onClick={() => this.applyRepresentation('default')} disabled={this.state.busy} >
                        Default
                    </Button>
                    <Button title='Applies cartoon polymer and ball-and-stick ligand representation preset'
                        onClick={() => this.applyRepresentation('cartoon')} disabled={this.state.busy} >
                        Cartoon
                    </Button>
                    <Button title='Applies spacefill representation preset'
                        onClick={() => this.applyRepresentation('spacefill')} disabled={this.state.busy} >
                        Spacefill
                    </Button>
                    <Button title='Applies molecular surface representation preset'
                        onClick={() => this.applyRepresentation('surface')} disabled={this.state.busy} >
                        Surface
                    </Button>
                </div>
            </NoncollapsableGroup>
            <NoncollapsableGroup title='Apply Style'>
                <div className='msp-flex-row'>
                    <Button title='Applies default appearance (no outline, no ignore-light)'
                        onClick={() => this.applyStyle('default')} disabled={this.state.busy} >
                        Default
                    </Button>
                    <Button title='Applies illustrative appearance (outline, ignore-light)'
                        onClick={() => this.applyStyle('illustrative')} disabled={this.state.busy} >
                        Illustrative
                    </Button>
                </div>
            </NoncollapsableGroup>
        </>;
    }
}


/** Visually imitates `ControlGroup` but is always expanded */
function NoncollapsableGroup(props: { title: string, children: any }): JSX.Element {
    return <div className='msp-control-group-wrapper'>
        <div className='msp-control-group-header'><div><b>{props.title}</b></div></div>
        {props.children}
    </div>;
}

async function applyRepresentationPreset(plugin: PluginContext, preset: PresetName) {
    const { structures } = plugin.managers.structure.hierarchy.selection;

    switch (preset) {
        case 'default':
            const defaultPreset = plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
            const provider = plugin.builders.structure.representation.resolveProvider(defaultPreset);
            await plugin.managers.structure.component.applyPreset(structures, provider);
            break;
        case 'spacefill':
            await plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations.illustrative);
            break;
        case 'cartoon':
            await plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations['polymer-and-ligand']);
            break;
        case 'surface':
            await plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations['molecular-surface']);
            break;
    }
}

async function applyStyle(plugin: PluginContext, style: StyleName) {
    if (style === 'default') {
        await plugin.managers.structure.component.setOptions({ ...plugin.managers.structure.component.state.options, ignoreLight: false });

        if (plugin.canvas3d) {
            const p = PD.getDefaultValues(PostprocessingParams);
            plugin.canvas3d.setProps({
                postprocessing: { outline: p.outline, occlusion: p.occlusion, shadow: p.shadow }
            });
        }
    }

    if (style === 'illustrative') {
        await plugin.managers.structure.component.setOptions({ ...plugin.managers.structure.component.state.options, ignoreLight: true });

        if (plugin.canvas3d) {
            const pp = plugin.canvas3d.props.postprocessing;
            plugin.canvas3d.setProps({
                postprocessing: {
                    outline: {
                        name: 'on',
                        params: pp.outline.name === 'on'
                            ? pp.outline.params
                            : {
                                scale: 1,
                                color: Color(0x000000),
                                threshold: 0.33,
                                includeTransparent: true,
                            }
                    },
                    occlusion: {
                        name: 'on',
                        params: pp.occlusion.name === 'on'
                            ? pp.occlusion.params
                            : {
                                multiScale: { name: 'off', params: {} },
                                radius: 5,
                                bias: 0.8,
                                blurKernelSize: 15,
                                blurDepthBias: 0.5,
                                samples: 32,
                                resolutionScale: 1,
                                color: Color(0x000000),
                                transparentThreshold: 0.4,
                            }
                    },
                    shadow: { name: 'off', params: {} },
                }
            });
        }
    }
}
