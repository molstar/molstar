/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AssemblySymmetryConfig } from '../../extensions/assembly-symmetry';
import { DefaultPluginUISpec, PluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginBehaviors } from '../../mol-plugin/behavior';
import { PluginConfig } from '../../mol-plugin/config';
import { ExtensionMap } from './extensions';
import { DefaultViewerOptions, ViewerOptions } from './options';
import { NoPrimaryFocusLociBindings } from '../../mol-plugin/behavior/dynamic/camera';
import { PluginSpec } from '../../mol-plugin/spec';
import { ViewerAutoPreset } from './presets';

export function createViewerSpec(options: Partial<ViewerOptions> = {}): PluginUISpec {
    const definedOptions = {} as any;
    // filter for defined properies only so the default values
    // are property applied
    for (const p of Object.keys(options) as (keyof ViewerOptions)[]) {
        if (options[p] !== void 0) definedOptions[p] = options[p];
    }

    const o: ViewerOptions = { ...DefaultViewerOptions, ...definedOptions };
    const defaultSpec = DefaultPluginUISpec();

    const disabledExtension = new Set(o.disabledExtensions ?? []);
    let baseBehaviors = defaultSpec.behaviors;

    if (o.viewportFocusBehavior === 'disabled') {
        baseBehaviors = baseBehaviors.filter(b =>
            b.transformer !== PluginBehaviors.Camera.FocusLoci
            && b.transformer !== PluginBehaviors.Representation.FocusLoci
        );
    } else if (o.viewportFocusBehavior === 'secondary-zoom') {
        baseBehaviors = baseBehaviors.filter(b =>
            b.transformer !== PluginBehaviors.Camera.FocusLoci
            && b.transformer !== PluginBehaviors.Representation.FocusLoci
        );

        baseBehaviors.push(PluginSpec.Behavior(PluginBehaviors.Camera.FocusLoci, {
            bindings: NoPrimaryFocusLociBindings
        }));
    }

    const spec: PluginUISpec = {
        canvas3d: {
            ...defaultSpec.canvas3d,
        },
        actions: defaultSpec.actions,
        behaviors: [
            ...baseBehaviors,
            ...o.extensions.filter(e => !disabledExtension.has(e)).map(e => ExtensionMap[e]),
        ],
        animations: [...defaultSpec.animations || []],
        customParamEditors: defaultSpec.customParamEditors,
        customFormats: o?.customFormats,
        layout: {
            initial: {
                isExpanded: o.layoutIsExpanded,
                showControls: o.layoutShowControls,
                controlsDisplay: o.layoutControlsDisplay,
                regionState: {
                    bottom: 'full',
                    left: o.collapseLeftPanel ? 'collapsed' : 'full',
                    right: o.collapseRightPanel ? 'hidden' : 'full',
                    top: 'full',
                }
            },
        },
        components: {
            ...defaultSpec.components,
            controls: {
                ...defaultSpec.components?.controls,
                top: o.layoutShowSequence ? undefined : 'none',
                bottom: o.layoutShowLog ? undefined : 'none',
                left: o.layoutShowLeftPanel ? undefined : 'none',
            },
            remoteState: o.layoutShowRemoteState ? 'default' : 'none',
        },
        config: [
            [PluginConfig.General.DisableAntialiasing, o.disableAntialiasing],
            [PluginConfig.General.PixelScale, o.pixelScale],
            [PluginConfig.General.PickScale, o.pickScale],
            [PluginConfig.General.Transparency, o.transparency],
            [PluginConfig.General.PreferWebGl1, o.preferWebgl1],
            [PluginConfig.General.AllowMajorPerformanceCaveat, o.allowMajorPerformanceCaveat],
            [PluginConfig.General.PowerPreference, o.powerPreference],
            [PluginConfig.General.ResolutionMode, o.resolutionMode],
            [PluginConfig.Viewport.ShowReset, o.viewportShowReset],
            [PluginConfig.Viewport.ShowScreenshotControls, o.viewportShowScreenshotControls],
            [PluginConfig.Viewport.ShowControls, o.viewportShowControls],
            [PluginConfig.Viewport.ShowExpand, o.viewportShowExpand],
            [PluginConfig.Viewport.ShowToggleFullscreen, o.viewportShowToggleFullscreen],
            [PluginConfig.Viewport.ShowSettings, o.viewportShowSettings],
            [PluginConfig.Viewport.ShowSelectionMode, o.viewportShowSelectionMode],
            [PluginConfig.Viewport.ShowAnimation, o.viewportShowAnimation],
            [PluginConfig.Viewport.ShowTrajectoryControls, o.viewportShowTrajectoryControls],
            [PluginConfig.State.DefaultServer, o.pluginStateServer],
            [PluginConfig.State.CurrentServer, o.pluginStateServer],
            [PluginConfig.VolumeStreaming.DefaultServer, o.volumeStreamingServer],
            [PluginConfig.VolumeStreaming.Enabled, !o.volumeStreamingDisabled],
            [PluginConfig.Download.DefaultPdbProvider, o.pdbProvider],
            [PluginConfig.Download.DefaultEmdbProvider, o.emdbProvider],
            [PluginConfig.Structure.DefaultRepresentationPreset, ViewerAutoPreset.id],
            [PluginConfig.Structure.SaccharideCompIdMapType, o.saccharideCompIdMapType],
            [AssemblySymmetryConfig.DefaultServerType, o.rcsbAssemblySymmetryDefaultServerType],
            [AssemblySymmetryConfig.DefaultServerUrl, o.rcsbAssemblySymmetryDefaultServerUrl],
            [AssemblySymmetryConfig.ApplyColors, o.rcsbAssemblySymmetryApplyColors],
            ...(o.config ?? []),
        ]
    };

    return spec;
}