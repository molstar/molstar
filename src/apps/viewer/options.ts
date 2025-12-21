/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssemblySymmetryConfig } from '../../extensions/assembly-symmetry';
import { G3dProvider } from '../../extensions/g3d/format';
import { SaccharideCompIdMapType } from '../../mol-model/structure/structure/carbohydrates/constants';
import { DataFormatProvider } from '../../mol-plugin-state/formats/provider';
import { PluginConfig, PluginConfigItem } from '../../mol-plugin/config';
import { PluginLayoutControlsDisplay } from '../../mol-plugin/layout';
import '../../mol-util/polyfill';
import { ObjectKeys } from '../../mol-util/type-helpers';
import { ExtensionMap } from './extensions';

const CustomFormats: [string, DataFormatProvider][] = [
    ['g3d', G3dProvider] as const
];

export const DefaultViewerOptions = {
    customFormats: CustomFormats as [string, DataFormatProvider][],
    extensions: ObjectKeys(ExtensionMap),
    disabledExtensions: [] as string[],
    layoutIsExpanded: true,
    layoutShowControls: true,
    layoutShowRemoteState: true,
    layoutControlsDisplay: 'reactive' as PluginLayoutControlsDisplay,
    layoutShowSequence: true,
    layoutShowLog: true,
    layoutShowLeftPanel: true,
    collapseLeftPanel: false,
    collapseRightPanel: false,
    disableAntialiasing: PluginConfig.General.DisableAntialiasing.defaultValue,
    pixelScale: PluginConfig.General.PixelScale.defaultValue,
    pickScale: PluginConfig.General.PickScale.defaultValue,
    transparency: PluginConfig.General.Transparency.defaultValue,
    preferWebgl1: PluginConfig.General.PreferWebGl1.defaultValue,
    allowMajorPerformanceCaveat: PluginConfig.General.AllowMajorPerformanceCaveat.defaultValue,
    powerPreference: PluginConfig.General.PowerPreference.defaultValue,
    resolutionMode: PluginConfig.General.ResolutionMode.defaultValue,
    illumination: false,

    viewportShowReset: PluginConfig.Viewport.ShowReset.defaultValue,
    viewportShowScreenshotControls: PluginConfig.Viewport.ShowScreenshotControls.defaultValue,
    viewportShowControls: PluginConfig.Viewport.ShowControls.defaultValue,
    viewportShowExpand: PluginConfig.Viewport.ShowExpand.defaultValue,
    viewportShowToggleFullscreen: PluginConfig.Viewport.ShowToggleFullscreen.defaultValue,
    viewportShowSettings: PluginConfig.Viewport.ShowSettings.defaultValue,
    viewportShowSelectionMode: PluginConfig.Viewport.ShowSelectionMode.defaultValue,
    viewportShowAnimation: PluginConfig.Viewport.ShowAnimation.defaultValue,
    viewportShowTrajectoryControls: PluginConfig.Viewport.ShowTrajectoryControls.defaultValue,
    viewportFocusBehavior: 'default' as 'default' | 'disabled',
    viewportBackgroundColor: undefined as string | undefined,

    pluginStateServer: PluginConfig.State.DefaultServer.defaultValue,
    volumeStreamingServer: PluginConfig.VolumeStreaming.DefaultServer.defaultValue,
    volumeStreamingDisabled: !PluginConfig.VolumeStreaming.Enabled.defaultValue,
    pdbProvider: PluginConfig.Download.DefaultPdbProvider.defaultValue,
    emdbProvider: PluginConfig.Download.DefaultEmdbProvider.defaultValue,
    saccharideCompIdMapType: 'default' as SaccharideCompIdMapType,
    rcsbAssemblySymmetryDefaultServerType: AssemblySymmetryConfig.DefaultServerType.defaultValue,
    rcsbAssemblySymmetryDefaultServerUrl: AssemblySymmetryConfig.DefaultServerUrl.defaultValue,
    rcsbAssemblySymmetryApplyColors: AssemblySymmetryConfig.ApplyColors.defaultValue,

    config: [] as [PluginConfigItem, any][],
};
export type ViewerOptions = typeof DefaultViewerOptions;