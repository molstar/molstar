/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PDBeStructureQualityReport } from '../../extensions/pdbe';
import { EmptyLoci } from '../../mol-model/loci';
import { StructureSelection } from '../../mol-model/structure';
import { AnimateModelIndex } from '../../mol-plugin-state/animation/built-in/model-index';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { createPluginUI } from '../../mol-plugin-ui';
import { PluginUIContext } from '../../mol-plugin-ui/context';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginCommands } from '../../mol-plugin/commands';
import { Script } from '../../mol-script/script';
import { Asset } from '../../mol-util/assets';
import { Color } from '../../mol-util/color';
import { StripedResidues } from './coloring';
import { CustomToastMessage } from './controls';
import { CustomColorThemeProvider } from './custom-theme';
import './index.html';
import { buildStaticSuperposition, dynamicSuperpositionTest, StaticSuperpositionTestData } from './superposition';
import '../../mol-plugin-ui/skin/light.scss';

type LoadParams = { url: string, format?: BuiltInTrajectoryFormat, isBinary?: boolean, assemblyId?: string }

class BasicWrapper {
    plugin: PluginUIContext;

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
                    }
                },
                components: {
                    remoteState: 'none'
                }
            }
        });

        this.plugin.representation.structure.themes.colorThemeRegistry.add(StripedResidues.colorThemeProvider!);
        this.plugin.representation.structure.themes.colorThemeRegistry.add(CustomColorThemeProvider);
        this.plugin.managers.lociLabels.addProvider(StripedResidues.labelProvider!);
        this.plugin.customModelProperties.register(StripedResidues.propertyProvider, true);

        this.plugin.managers.dragAndDrop.addHandler('custom-wrapper', (files) => {
            if (files.some(f => f.name.toLowerCase().endsWith('.testext'))) {
                console.log('.testext File dropped');
                return true;
            }
            return false;
        });
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
                params: {}
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

        const trackball = this.plugin.canvas3d.props.trackball;
        PluginCommands.Canvas3D.SetSettings(this.plugin, {
            settings: {
                trackball: {
                    ...trackball,
                    animate: trackball.animate.name === 'spin'
                        ? { name: 'off', params: {} }
                        : { name: 'spin', params: { speed: 1 } }
                }
            }
        });
        if (this.plugin.canvas3d.props.trackball.animate.name !== 'spin') {
            PluginCommands.Camera.Reset(this.plugin, {});
        }
    }

    private animateModelIndexTargetFps() {
        return Math.max(1, this.animate.modelIndex.targetFps | 0);
    }

    animate = {
        modelIndex: {
            targetFps: 8,
            onceForward: () => { this.plugin.managers.animation.play(AnimateModelIndex, { duration: { name: 'computed', params: { targetFps: this.animateModelIndexTargetFps() } }, mode: { name: 'once', params: { direction: 'forward' } } }); },
            onceBackward: () => { this.plugin.managers.animation.play(AnimateModelIndex, { duration: { name: 'computed', params: { targetFps: this.animateModelIndexTargetFps() } }, mode: { name: 'once', params: { direction: 'backward' } } }); },
            palindrome: () => { this.plugin.managers.animation.play(AnimateModelIndex, { duration: { name: 'computed', params: { targetFps: this.animateModelIndexTargetFps() } }, mode: { name: 'palindrome', params: {} } }); },
            loop: () => { this.plugin.managers.animation.play(AnimateModelIndex, { duration: { name: 'computed', params: { targetFps: this.animateModelIndexTargetFps() } }, mode: { name: 'loop', params: { direction: 'forward' } } }); },
            stop: () => this.plugin.managers.animation.stop()
        }
    };

    coloring = {
        applyStripes: async () => {
            this.plugin.dataTransaction(async () => {
                for (const s of this.plugin.managers.structure.hierarchy.current.structures) {
                    await this.plugin.managers.structure.component.updateRepresentationsTheme(s.components, { color: StripedResidues.propertyProvider.descriptor.name as any });
                }
            });
        },
        applyCustomTheme: async () => {
            this.plugin.dataTransaction(async () => {
                for (const s of this.plugin.managers.structure.hierarchy.current.structures) {
                    await this.plugin.managers.structure.component.updateRepresentationsTheme(s.components, { color: CustomColorThemeProvider.name as any });
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
    };

    interactivity = {
        highlightOn: () => {
            const data = this.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
            if (!data) return;

            const seq_id = 7;
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
    };

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
    };
}

(window as any).BasicMolStarWrapper = new BasicWrapper();