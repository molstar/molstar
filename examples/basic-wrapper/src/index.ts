/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PDBeStructureQualityReport } from '@molstar/pdbe-extension';
import { EmptyLoci } from '@molstar/model/model/loci';
import { StructureSelection } from '@molstar/model/model/structure';
import { AnimateModelIndex } from '@molstar/plugin/state/animation/built-in/model-index';
import type { BuiltInTrajectoryFormat } from '@molstar/plugin/state/formats/trajectory';
import { createPluginUI } from '@molstar/plugin-ui';
import { PluginUIContext } from '@molstar/plugin-ui/context';
import { renderReact18 } from '@molstar/plugin-ui/react18';
import { DefaultPluginUISpec } from '@molstar/plugin-ui/spec';
import { PluginCommands } from '@molstar/plugin/commands';
import { Script } from '@molstar/model/script/script';
import { Asset } from '@molstar/core/util/assets';
import { Color } from '@molstar/core/util/color';
import { StripedResidues } from '@molstar/basic-wrapper-example/coloring';
import { CustomToastMessage } from '@molstar/basic-wrapper-example/controls';
import { CustomColorThemeProvider } from '@molstar/basic-wrapper-example/custom-theme';
import './index.html';
import './tm-align.html';
import { buildStaticSuperposition, dynamicSuperpositionTest, StaticSuperpositionTestData, tmAlignStructures, loadStructuresNoAlignment, sequenceAlignStructures } from '@molstar/basic-wrapper-example/superposition';
import '@molstar/plugin-ui/skin/light.scss';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';

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
                        : { name: 'spin', params: { speed: 0.1, axis: Vec3.create(0, -1, 0) } }
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

    /**
     * Run TM-align on two structures
     * @param pdbId1 - PDB ID of first structure (reference)
     * @param chain1 - Chain ID of first structure
     * @param pdbId2 - PDB ID of second structure (mobile)
     * @param chain2 - Chain ID of second structure
     * @param color1 - Optional color for first structure (hex, default blue)
     * @param color2 - Optional color for second structure (hex, default red)
     */
    tmAlign(pdbId1: string, chain1: string, pdbId2: string, chain2: string, color1?: number, color2?: number) {
        return tmAlignStructures(this.plugin, pdbId1, chain1, pdbId2, chain2, color1, color2);
    }

    /**
     * Load two structures without alignment
     * @param pdbId1 - PDB ID of first structure
     * @param chain1 - Chain ID of first structure
     * @param pdbId2 - PDB ID of second structure
     * @param chain2 - Chain ID of second structure
     * @param color1 - Optional color for first structure (hex, default blue)
     * @param color2 - Optional color for second structure (hex, default red)
     */
    loadStructures(pdbId1: string, chain1: string, pdbId2: string, chain2: string, color1?: number, color2?: number) {
        return loadStructuresNoAlignment(this.plugin, pdbId1, chain1, pdbId2, chain2, color1, color2);
    }

    /**
     * Align two structures using sequence alignment
     * @param pdbId1 - PDB ID of first structure (reference)
     * @param chain1 - Chain ID of first structure
     * @param pdbId2 - PDB ID of second structure (mobile)
     * @param chain2 - Chain ID of second structure
     * @param color1 - Optional color for first structure (hex, default blue)
     * @param color2 - Optional color for second structure (hex, default red)
     */
    sequenceAlign(pdbId1: string, chain1: string, pdbId2: string, chain2: string, color1?: number, color2?: number) {
        return sequenceAlignStructures(this.plugin, pdbId1, chain1, pdbId2, chain2, color1, color2);
    }
}

(window as any).BasicMolStarWrapper = new BasicWrapper();