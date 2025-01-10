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
import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { CollapsableControls, PurePluginUIComponent } from '../base';
import { Button, IconButton } from '../controls/common';
import { Draw, MagicWandSvg } from '../controls/icons';

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

interface QuickStylesState {
    preset: 'default' | 'spacefill' | 'surface',
    stylized: boolean,
}

export class QuickStyles extends PurePluginUIComponent<{}, QuickStylesState> {
    state: QuickStylesState = { preset: 'default', stylized: false };

    async default() {
        const { structures } = this.plugin.managers.structure.hierarchy.selection;
        const preset = this.plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
        const provider = this.plugin.builders.structure.representation.resolveProvider(preset);
        await this.plugin.managers.structure.component.applyPreset(structures, provider);
        await this.setStylized(this.state.stylized);
        this.setState({ preset: 'default' });
    }

    async illustrative() {
        const { structures } = this.plugin.managers.structure.hierarchy.selection;
        await this.plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations.illustrative);
        await this.setStylized(this.state.stylized);
        this.setState({ preset: 'spacefill' });
    }
    async surface() {
        const { structures } = this.plugin.managers.structure.hierarchy.selection;
        await this.plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations['molecular-surface']);
        await this.setStylized(this.state.stylized);
        this.setState({ preset: 'surface' });
    }

    async setStylized(stylized: boolean) {
        if (stylized) {
            this.plugin.managers.structure.component.setOptions({ ...this.plugin.managers.structure.component.state.options, ignoreLight: true });

            if (this.plugin.canvas3d) {
                const pp = this.plugin.canvas3d.props.postprocessing;
                this.plugin.canvas3d.setProps({
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
        } else {
            this.plugin.managers.structure.component.setOptions({ ...this.plugin.managers.structure.component.state.options, ignoreLight: false });

            if (this.plugin.canvas3d) {
                const p = PD.getDefaultValues(PostprocessingParams);
                this.plugin.canvas3d.setProps({
                    postprocessing: { outline: p.outline, occlusion: p.occlusion }
                });
            }
        }
        this.setState({ stylized });
    }

    async toggleStylized() {
        await this.setStylized(!this.state.stylized);
    }

    render() {
        return <div className='msp-flex-row'>
            <Button noOverflow title='Applies default representation preset' onClick={() => this.default()} style={{ width: 'auto' }}>
                Default
            </Button>
            <Button noOverflow title='Applies illustrative representation preset' onClick={() => this.illustrative()} style={{ width: 'auto' }}>
                Spacefill
            </Button>
            <Button noOverflow title='Applies molecular surface representation preset' onClick={() => this.surface()} style={{ width: 'auto' }}>
                Surface
            </Button>
            <IconButton
                title='Stylize. &#10;Does not change representation, toggles stylized appearance (outline, occlusion effects, ignore-light representation parameter).'
                onClick={() => this.toggleStylized()} style={{ width: 'auto', marginInline: 10 }} svg={Draw} toggleState={this.state.stylized}
            ></IconButton>
        </div>;
    }
}
