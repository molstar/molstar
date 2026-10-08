/**
 * Copyright (c) 2022-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mp4Export } from '@molstar/mp4-export-extension';
import { Mmcif } from '@molstar/plugin/state/formats/trajectory/mmcif';
import { Spacefill } from '@molstar/plugin/registry/structure/spacefill';
import { UniformColorThemeProvider } from '@molstar/graphics/theme/color/uniform';
import { IllustrativeColorThemeProvider } from '@molstar/graphics/theme/color/illustrative';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { createPluginUI } from '@molstar/plugin-ui';
import { renderReact18 } from '@molstar/plugin-ui/react18';
import { PluginUIContext } from '@molstar/plugin-ui/context';
import { DefaultPluginUIComponents, DefaultPluginUICustomParamEditors } from '@molstar/plugin-ui/default-ui';
import type { PluginUISpec } from '@molstar/plugin-ui/spec';
import { PluginConfig } from '@molstar/plugin/config';
import type { PluginLayoutControlsDisplay } from '@molstar/plugin/layout';
import { PluginSpec, type PluginRegistryEntry } from '@molstar/plugin/spec';
import '@molstar/core/util/polyfill';
import { ObjectKeys } from '@molstar/core/util/type-helpers';
import type { SaccharideCompIdMapType } from '@molstar/model/model/structure/structure/carbohydrates/constants';
import { Backgrounds } from '@molstar/backgrounds-extension';
import { LeftPanel, RightPanel } from '@molstar/mesoscale-explorer/ui/panels';
import { Color } from '@molstar/core/util/color';
import { PluginBehaviors } from '@molstar/plugin/behavior';
import { MesoFocusLoci } from '@molstar/mesoscale-explorer/behavior/camera';
import { type GraphicsMode, MesoscaleState } from '@molstar/mesoscale-explorer/data/state';
import { MesoSelectLoci } from '@molstar/mesoscale-explorer/behavior/select';
import type { Transparency } from '@molstar/graphics/gl/webgl/render-item';
import {
  LoadModel,
  loadExampleEntry,
  loadPdb,
  loadPdbIhm,
  loadUrl,
  openState,
} from '@molstar/mesoscale-explorer/ui/states';
import { Asset } from '@molstar/core/util/assets';
import { AnimateCameraSpin } from '@molstar/plugin/state/animation/built-in/camera-spin';
import { AnimateCameraRock } from '@molstar/plugin/state/animation/built-in/camera-rock';
import { AnimateStateSnapshots } from '@molstar/plugin/state/animation/built-in/state-snapshots';
import { MesoViewportSnapshotDescription } from '@molstar/mesoscale-explorer/ui/entities';

export { PLUGIN_VERSION as version } from '@molstar/plugin/version';
export { setDebugMode, setProductionMode, setTimingMode, consoleStats } from '@molstar/core/util/debug';

export type ExampleEntry = {
  id: string;
  label: string;
  url: string;
  type: 'molx' | 'molj' | 'cif' | 'bcif';
  description?: string;
  link?: string;
};

export type MesoscaleExplorerState = {
  examples?: ExampleEntry[];
  graphicsMode: GraphicsMode;
  illumination: boolean;
  stateRef?: string;
  driver?: any;
  stateCache: { [k: string]: any };
};

//

/**
 * The providers the explorer lists, and nothing else. Its presets set themes through direct state updates, so every
 * theme it names must be here, or the registry default is substituted silently: the spacefill representation with its
 * element-symbol color and `physical` size themes (the explorer names `physical`), and the `uniform` and `illustrative`
 * structure color themes. It lists the mmCIF format with its actions because it reads `.cif` and `.bcif` and needs
 * their extensions registered (the other formats of a zip archive are parsed through transformers directly), its three
 * animations, and the custom formats, which come last.
 *
 * As in 5.x, a custom name that equals a built-in format name overrides it: the mmCIF entry is replaced by a copy
 * without that provider, so no name is registered twice.
 */
function createRegistry(customFormats: [string, DataFormatProvider.Unnamed][] | undefined): PluginRegistryEntry[] {
  const formats = (customFormats ?? []).map(([name, provider]) => DataFormatProvider.withName(provider, name));
  const customNames = new Set(formats.map((f) => f.name));
  return [
    Spacefill,
    { structure: { themes: { color: [UniformColorThemeProvider, IllustrativeColorThemeProvider] } } },
    { ...Mmcif, formats: Mmcif.formats!.filter((f) => !customNames.has(f.name)) },
    { animations: [AnimateCameraSpin, AnimateCameraRock, AnimateStateSnapshots] },
    { formats },
  ];
}

const Extensions = {
  backgrounds: PluginSpec.Behavior(Backgrounds),
  'mp4-export': PluginSpec.Behavior(Mp4Export),
};

const DefaultMesoscaleExplorerOptions = {
  customFormats: [] as [string, DataFormatProvider.Unnamed][],
  extensions: ObjectKeys(Extensions),
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
  transparency: 'blended' as Transparency,
  preferWebgl1: PluginConfig.General.PreferWebGl1.defaultValue,
  allowMajorPerformanceCaveat: PluginConfig.General.AllowMajorPerformanceCaveat.defaultValue,
  powerPreference: PluginConfig.General.PowerPreference.defaultValue,
  resolutionMode: PluginConfig.General.ResolutionMode.defaultValue,
  illumination: false,

  viewportShowExpand: PluginConfig.Viewport.ShowExpand.defaultValue,
  viewportShowControls: PluginConfig.Viewport.ShowControls.defaultValue,
  viewportShowSettings: PluginConfig.Viewport.ShowSettings.defaultValue,
  viewportShowSelectionMode: true,
  viewportShowAnimation: false,
  viewportShowTrajectoryControls: false,
  pluginStateServer: PluginConfig.State.DefaultServer.defaultValue,
  volumeStreamingServer: PluginConfig.VolumeStreaming.DefaultServer.defaultValue,
  volumeStreamingDisabled: !PluginConfig.VolumeStreaming.Enabled.defaultValue,
  pdbProvider: PluginConfig.Download.DefaultPdbProvider.defaultValue,
  emdbProvider: PluginConfig.Download.DefaultEmdbProvider.defaultValue,
  saccharideCompIdMapType: 'default' as SaccharideCompIdMapType,

  graphicsMode: 'quality' as GraphicsMode,
  driver: undefined,
};
type MesoscaleExplorerOptions = typeof DefaultMesoscaleExplorerOptions;

export class MesoscaleExplorer {
  constructor(public plugin: PluginUIContext) {}

  async loadExample(id: string) {
    const entries = (this.plugin.customState as MesoscaleExplorerState).examples || [];
    const entry = entries.find((e) => e.id === id);
    if (entry !== undefined) {
      await loadExampleEntry(this.plugin, entry);
    }
  }

  async loadUrl(url: string, type: 'molx' | 'molj' | 'cif' | 'bcif') {
    await loadUrl(this.plugin, url, type);
  }

  async loadPdb(id: string) {
    await loadPdb(this.plugin, id);
  }

  /**
   * @deprecated Scheduled for removal in v5. Use {@link loadPdbIhm | loadPdbIhm(id: string)} instead.
   */
  async loadPdbDev(id: string) {
    await this.loadPdbIhm(id);
  }

  async loadPdbIhm(id: string) {
    await loadPdbIhm(this.plugin, id);
  }

  static async create(elementOrId: string | HTMLElement, options: Partial<MesoscaleExplorerOptions> = {}) {
    const definedOptions = {} as any;
    // filter for defined properies only so the default values
    // are property applied
    for (const p of Object.keys(options) as (keyof MesoscaleExplorerOptions)[]) {
      if (options[p] !== void 0) definedOptions[p] = options[p];
    }

    const o: MesoscaleExplorerOptions = { ...DefaultMesoscaleExplorerOptions, ...definedOptions };
    const defaultComponents = DefaultPluginUIComponents();

    const spec: PluginUISpec = {
      registry: createRegistry(o.customFormats),
      behaviors: [
        PluginSpec.Behavior(PluginBehaviors.Camera.CameraAxisHelper),
        PluginSpec.Behavior(PluginBehaviors.Camera.CameraControls),
        PluginSpec.Behavior(PluginBehaviors.State.SnapshotControls),

        PluginSpec.Behavior(MesoFocusLoci),
        PluginSpec.Behavior(MesoSelectLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.SelectLoci),
        ...o.extensions.map((e) => Extensions[e]),
      ],
      customParamEditors: DefaultPluginUICustomParamEditors(),
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
          },
        },
      },
      components: {
        ...defaultComponents,
        controls: {
          ...defaultComponents.controls,
          top: 'none',
          bottom: 'none',
          left: LeftPanel,
          right: RightPanel,
        },
        remoteState: 'none',
        viewport: {
          snapshotDescription: MesoViewportSnapshotDescription,
        },
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
        [PluginConfig.Viewport.ShowExpand, o.viewportShowExpand],
        [PluginConfig.Viewport.ShowControls, o.viewportShowControls],
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
        [PluginConfig.Structure.SaccharideCompIdMapType, o.saccharideCompIdMapType],
      ],
    };

    const element = typeof elementOrId === 'string' ? document.getElementById(elementOrId) : elementOrId;
    if (!element) throw new Error(`Could not get element with id '${elementOrId}'`);

    const plugin = await createPluginUI({
      target: element,
      spec,
      render: renderReact18,
      onBeforeUIRender: async (plugin) => {
        let examples: MesoscaleExplorerState['examples'] = undefined;
        try {
          examples = await plugin.fetch({ url: '../examples/list.json', type: 'json' }).run();
          // extend the array with file tour.json if it exists
          const tour = await plugin.fetch({ url: '../examples/tour.json', type: 'json' }).run();
          if (tour) {
            examples = examples?.concat(tour);
          }
        } catch (e) {
          console.log(e);
        }

        (plugin.customState as MesoscaleExplorerState) = {
          examples,
          graphicsMode: o.graphicsMode,
          illumination: o.illumination,
          driver: o.driver,
          stateCache: {},
        };

        await MesoscaleState.init(plugin);
      },
    });

    plugin.canvas3d?.setProps({
      renderer: {
        backgroundColor: Color(0x101010),
      },
      cameraFog: { name: 'off', params: {} },
      hiZ: { enabled: true },
      xr: {
        disablePostprocessing: false,
        sceneRadiusInMeters: 0.75,
      },
    });

    plugin.state.setSnapshotParams({
      image: true,
      componentManager: false,
      structureSelection: true,
    });

    plugin.managers.dragAndDrop.addHandler('mesoscale-explorer', (files) => {
      const sessions = files.filter((f) => {
        const fn = f.name.toLowerCase();
        return fn.endsWith('.molx') || fn.endsWith('.molj');
      });

      if (sessions.length > 0) {
        openState(plugin, sessions[0]);
      } else {
        plugin.runTask(
          plugin.state.data.applyAction(LoadModel, {
            files: files.map((f) => Asset.File(f)),
          }),
        );
      }

      return true;
    });

    plugin.state.events.object.created.subscribe((e) => {
      (plugin.customState as MesoscaleExplorerState).stateCache = {};
    });

    plugin.state.events.object.removed.subscribe((e) => {
      (plugin.customState as MesoscaleExplorerState).stateCache = {};
    });

    return new MesoscaleExplorer(plugin);
  }

  handleResize() {
    this.plugin.layout.events.updated.next(void 0);
  }

  dispose() {
    this.plugin.dispose();
  }
}
