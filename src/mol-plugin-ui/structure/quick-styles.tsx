/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PresetStructureRepresentations } from '../../mol-plugin-state/builder/structure/representation-preset';
import { Color } from '../../mol-util/color';
import { CollapsableControls, PurePluginUIComponent } from '../base';
import { Button } from '../controls/common';
import { MagicWandSvg } from '../controls/icons';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PostprocessingParams } from '../../mol-canvas3d/passes/postprocessing';
import { PluginConfig } from '../../mol-plugin/config';
import { StructureComponentManager } from '../../mol-plugin-state/manager/structure/component';

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

export class QuickStyles extends PurePluginUIComponent {
    async default() {
        const { structures } = this.plugin.managers.structure.hierarchy.selection;
        const preset = this.plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
        const provider = this.plugin.builders.structure.representation.resolveProvider(preset);
        await this.plugin.managers.structure.component.applyPreset(structures, provider);

        this.plugin.managers.structure.component.setOptions(PD.getDefaultValues(StructureComponentManager.OptionsParams));

        if (this.plugin.canvas3d) {
            const p = PD.getDefaultValues(PostprocessingParams);
            this.plugin.canvas3d.setProps({
                postprocessing: { outline: p.outline, occlusion: p.occlusion }
            });
        }
    }

    async illustrative() {
        const { structures } = this.plugin.managers.structure.hierarchy.selection;
        await this.plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations.illustrative);

        if (this.plugin.canvas3d) {
            this.plugin.canvas3d.setProps({
                postprocessing: {
                    outline: {
                        name: 'on',
                        params: { scale: 1, color: Color(0x000000), threshold: 0.25 }
                    },
                    occlusion: {
                        name: 'on',
                        params: { bias: 0.8, blurKernelSize: 15, radius: 5, samples: 32, resolutionScale: 1 }
                    },
                }
            });
        }
    }

    async stylized() {
        this.plugin.managers.structure.component.setOptions({ ...this.plugin.managers.structure.component.state.options, ignoreLight: true });

        if (this.plugin.canvas3d) {
            const pp = this.plugin.canvas3d.props.postprocessing;
            this.plugin.canvas3d.setProps({
                postprocessing: {
                    outline: {
                        name: 'on',
                        params: pp.outline.name === 'on'
                            ? pp.outline.params
                            : { scale: 1, color: Color(0x000000), threshold: 0.33 }
                    },
                    occlusion: {
                        name: 'on',
                        params: pp.occlusion.name === 'on'
                            ? pp.occlusion.params
                            : { bias: 0.8, blurKernelSize: 15, radius: 5, samples: 32, resolutionScale: 1 }
                    },
                }
            });
        }
    }

    render() {
        return <div className='msp-flex-row'>
            <Button noOverflow title='Applies default representation preset. Set outline and occlusion effects to defaults.' onClick={() => this.default()} style={{ width: 'auto' }}>
                Default
            </Button>
            <Button noOverflow title='Applies no representation preset. Enables outline and occlusion effects. Enables ignore-light representation parameter.' onClick={() => this.stylized()} style={{ width: 'auto' }}>
                Stylized
            </Button>
            <Button noOverflow title='Applies illustrative representation preset. Enables outline and occlusion effects. Enables ignore-light parameter.' onClick={() => this.illustrative()} style={{ width: 'auto' }}>
                Illustrative
            </Button>
        </div>;
    }
}