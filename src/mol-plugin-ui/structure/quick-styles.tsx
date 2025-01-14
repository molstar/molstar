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
import { ToggleButton } from '../controls/common';
import { MagicWandSvg } from '../controls/icons';
import { ParameterControls } from '../controls/parameters';


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


const QuickStyleParams = {
    preset: PD.Select('', [['', ''], ['default', 'Default'], ['cartoon', 'Cartoon'], ['spacefill', 'Spacefill'], ['surface', 'Surface']] as const, { isHidden: true }),
    stylized: PD.Boolean(false, { description: 'Toggles stylized appearance (outline, occlusion effects, ignore-light representation parameter). Does not affect representation type.' }),
};

type QuickStyleProps = PD.Values<typeof QuickStyleParams>;

interface QuickStylesState {
    props: QuickStyleProps,
    busy: boolean,
}


export class QuickStyles extends PurePluginUIComponent<{}, QuickStylesState> {
    state: QuickStylesState = { props: { preset: '', stylized: false }, busy: false };

    componentDidMount() {
        let presetInvalid = false;
        const syncPreset = () => {
            if (presetInvalid) {
                presetInvalid = false;
                this.setState(old => ({ props: { ...old.props, preset: '' } }));
            }
        };
        this.subscribe(this.plugin.state.data.events.changed, e => {
            presetInvalid = true;
            if (!e.inTransaction && !this.plugin.behaviors.state.isAnimating.value) {
                syncPreset();
            } // else: wait for animation to finish, then syncs preset
        });
        this.subscribe(this.plugin.behaviors.state.isAnimating, isAnimating => {
            if (!isAnimating && !this.plugin.behaviors.state.isUpdating.value) {
                syncPreset();
            }
        });
    }

    async setPreset(preset: QuickStyleProps['preset']) {
        const { structures } = this.plugin.managers.structure.hierarchy.selection;

        switch (preset) {
            case 'default':
                const defaultPreset = this.plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
                const provider = this.plugin.builders.structure.representation.resolveProvider(defaultPreset);
                await this.plugin.managers.structure.component.applyPreset(structures, provider);
                break;
            case 'spacefill':
                await this.plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations.illustrative);
                break;
            case 'cartoon':
                await this.plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations['polymer-and-ligand']);
                break;
            case 'surface':
                await this.plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations['molecular-surface']);
                break;
        }
        await this.setStylized(this.state.props.stylized);
    }

    async setStylized(stylized: boolean) {
        if (stylized) {
            await this.plugin.managers.structure.component.setOptions({ ...this.plugin.managers.structure.component.state.options, ignoreLight: true });

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
            await this.plugin.managers.structure.component.setOptions({ ...this.plugin.managers.structure.component.state.options, ignoreLight: false });

            if (this.plugin.canvas3d) {
                const p = PD.getDefaultValues(PostprocessingParams);
                this.plugin.canvas3d.setProps({
                    postprocessing: { outline: p.outline, occlusion: p.occlusion }
                });
            }
        }
    }

    async changeValues(props: QuickStyleProps) {
        const currentProps = this.state.props;
        this.setState({ busy: true });
        if (props.preset !== currentProps.preset) {
            await this.setPreset(props.preset);
        }
        if (props.stylized !== currentProps.stylized) {
            await this.setStylized(props.stylized);
        }
        this.setState({ props, busy: false });
    }

    async changePreset(preset: QuickStyleProps['preset']) {
        await this.changeValues({ ...this.state.props, preset });
    }

    render() {
        return <>
            <div className='msp-flex-row'>
                <ToggleButton label='Default' title='Applies default representation preset (depends on structure size)'
                    toggle={() => this.changePreset('default')} isSelected={this.state.props.preset === 'default'} disabled={this.state.busy} />
                <ToggleButton label='Cartoon' title='Applies cartoon polymer + ball-and-stick ligand representation preset'
                    toggle={() => this.changePreset('cartoon')} isSelected={this.state.props.preset === 'cartoon'} disabled={this.state.busy} />
                <ToggleButton label='Spacefill' title='Applies spacefill representation preset'
                    toggle={() => this.changePreset('spacefill')} isSelected={this.state.props.preset === 'spacefill'} disabled={this.state.busy} />
                <ToggleButton label='Surface' title='Applies molecular surface representation preset'
                    toggle={() => this.changePreset('surface')} isSelected={this.state.props.preset === 'surface'} disabled={this.state.busy} />
            </div>
            <ParameterControls params={QuickStyleParams} values={this.state.props} onChangeValues={props => this.changeValues(props)} isDisabled={this.state.busy} />
        </>;
    }
}
