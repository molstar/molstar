/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { EmptyLoci } from '../../mol-model/loci';
import { StructureSelection } from '../../mol-model/structure';
import { createPlugin, DefaultPluginSpec } from '../../mol-plugin';
import { AnimateModelIndex } from '../../mol-plugin-state/animation/built-in';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { Script } from '../../mol-script/script';
import { Color } from '../../mol-util/color';
import { StripedResidues } from './coloring';
import { CustomToastMessage } from './controls';
import './index.html';
import { buildStaticSuperposition, dynamicSuperpositionTest, StaticSuperpositionTestData } from './superposition';
import { PDBeStructureQualityReport } from '../../extensions/pdbe';
import { Asset } from '../../mol-util/assets';
require('mol-plugin-ui/skin/light.scss');

type LoadParams = { url: string, format?: BuiltInTrajectoryFormat, isBinary?: boolean, assemblyId?: string }

class BasicWrapper {
    plugin: PluginContext;

    init(target: string | HTMLElement) {
        this.plugin = createPlugin(typeof target === 'string' ? document.getElementById(target)! : target, {
            ...DefaultPluginSpec,
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: false
                },
                controls: {
                    // left: 'none'
                }
            },
            components: {
                remoteState: 'none'
            }
        });

        this.plugin.representation.structure.themes.colorThemeRegistry.add(StripedResidues.colorThemeProvider!);
        this.plugin.managers.lociLabels.addProvider(StripedResidues.labelProvider!);
        this.plugin.customModelProperties.register(StripedResidues.propertyProvider, true);
    }

    async load({ url, format = 'mmcif', isBinary = false, assemblyId = '' }: LoadParams) {
        await this.plugin.clear();

        const data = await this.plugin.builders.data.download({ url: Asset.Url(url), isBinary }, { state: { isGhost: true } });
        const trajectory = await this.plugin.builders.structure.parseTrajectory(data, format);

        await this.plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
            structure: assemblyId ? {
                name: 'assembly',
                params: { id: assemblyId }
            } : {
                name: 'model',
                params: { }
            },
            showUnitcell: false,
            representationPreset: 'auto'
        });
    }

    setBackground(color: number) {
        PluginCommands.Canvas3D.SetSettings(this.plugin, { settings: props => { props.renderer.backgroundColor = Color(color); } });
    }

    toggleSpin() {
        if (!this.plugin.canvas3d) return;

        PluginCommands.Canvas3D.SetSettings(this.plugin, {
            settings: props => {
                props.trackball.spin = !props.trackball.spin;
            }
        });
        if (!this.plugin.canvas3d.props.trackball.spin) PluginCommands.Camera.Reset(this.plugin, {});
    }

    animate = {
        modelIndex: {
            maxFPS: 8,
            onceForward: () => { this.plugin.managers.animation.play(AnimateModelIndex, { maxFPS: Math.max(0.5, this.animate.modelIndex.maxFPS | 0), mode: { name: 'once', params: { direction: 'forward' } } }); },
            onceBackward: () => { this.plugin.managers.animation.play(AnimateModelIndex, { maxFPS: Math.max(0.5, this.animate.modelIndex.maxFPS | 0), mode: { name: 'once', params: { direction: 'backward' } } }); },
            palindrome: () => { this.plugin.managers.animation.play(AnimateModelIndex, { maxFPS: Math.max(0.5, this.animate.modelIndex.maxFPS | 0), mode: { name: 'palindrome', params: {} } }); },
            loop: () => { this.plugin.managers.animation.play(AnimateModelIndex, { maxFPS: Math.max(0.5, this.animate.modelIndex.maxFPS | 0), mode: { name: 'loop', params: {} } }); },
            stop: () => this.plugin.managers.animation.stop()
        }
    }

    coloring = {
        applyStripes: async () => {
            this.plugin.dataTransaction(async () => {
                for (const s of this.plugin.managers.structure.hierarchy.current.structures) {
                    await this.plugin.managers.structure.component.updateRepresentationsTheme(s.components, { color: StripedResidues.propertyProvider.descriptor.name as any });
                }
            });
        },
        applyDefault: async () => {
            this.plugin.dataTransaction(async () => {
                for (const s of this.plugin.managers.structure.hierarchy.current.structures) {
                    await this.plugin.managers.structure.component.updateRepresentationsTheme(s.components, { color: 'default' });
                }
            });
        }
    }

    interactivity = {
        highlightOn: () => {
            const seq_id = 7;
            const data = (this.plugin.state.data.select('asm')[0].obj as PluginStateObject.Molecule.Structure).data;
            const sel = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
                'residue-test': Q.core.rel.eq([Q.struct.atomProperty.macromolecular.label_seq_id(), seq_id]),
                'group-by': Q.struct.atomProperty.macromolecular.residueKey()
            }), data);
            const loci = StructureSelection.toLociWithSourceUnits(sel);
            this.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci });
        },
        clearHighlight: () => {
            this.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci: EmptyLoci });
        }
    }

    tests = {
        staticSuperposition: async () => {
            await this.plugin.clear();
            return buildStaticSuperposition(this.plugin, StaticSuperpositionTestData);
        },
        dynamicSuperposition: async () => {
            await this.plugin.clear();
            return dynamicSuperpositionTest(this.plugin, ['1tqn', '2hhb', '4hhb'], 'HEM');
        },
        toggleValidationTooltip: () => {
            return this.plugin.state.updateBehavior(PDBeStructureQualityReport, params => { params.showTooltip = !params.showTooltip; });
        },
        showToasts: () => {
            PluginCommands.Toast.Show(this.plugin, {
                title: 'Toast 1',
                message: 'This is an example text, timeout 3s',
                key: 'toast-1',
                timeoutMs: 3000
            });
            PluginCommands.Toast.Show(this.plugin, {
                title: 'Toast 2',
                message: CustomToastMessage,
                key: 'toast-2'
            });
        },
        hideToasts: () => {
            PluginCommands.Toast.Hide(this.plugin, { key: 'toast-1' });
            PluginCommands.Toast.Hide(this.plugin, { key: 'toast-2' });
        }
    }
}

(window as any).BasicMolStarWrapper = new BasicWrapper();