/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import type { PluginStateAnimation } from '@molstar/plugin/state/animation/model';
import type { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import type { StateAction, StateTransformer } from '@molstar/core/state';
import type { PluginConfigItem } from '@molstar/plugin/config';
import type { PluginLayoutStateProps } from '@molstar/plugin/layout';
import type { ParticleRepresentationProvider } from '@molstar/graphics/repr/particles/representation';
import type { StructureRepresentationProvider } from '@molstar/graphics/repr/structure/representation';
import type { VolumeRepresentationProvider } from '@molstar/graphics/repr/volume/representation';
import type { ColorTheme } from '@molstar/graphics/theme/color';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import type { StructureSelectionQuery } from '@molstar/plugin/state/queries/structure/query';
import type { TrajectoryHierarchyPresetProvider } from '@molstar/plugin/state/builder/structure/hierarchy-presets/types';
import type { StructureRepresentationPresetProvider } from '@molstar/plugin/state/builder/structure/representation-presets/types';
import type { LociLabelProvider } from '@molstar/plugin/state/manager/loci-label';
import type { MarkdownExtension } from '@molstar/plugin/state/manager/markdown-extensions';
import type { PluginDragAndDropEntry } from '@molstar/plugin/state/manager/drag-and-drop';

export { PluginSpec };

interface PluginSpec {
  /** Providers registered by `PluginContext.init()`, in order, before behaviors are initialized. */
  registry?: readonly PluginRegistryEntry[];
  behaviors: PluginSpec.Behavior[];
  canvas3d?: PartialCanvas3DProps;
  layout?: {
    initial?: Partial<PluginLayoutStateProps>;
  };
  config?: [PluginConfigItem, unknown][];
}

/**
 * A plain record of providers grouped by registry. Entries are pure data: providers only, with no behaviors, config,
 * or functions run at registration. See `PluginContext.register`.
 */
export interface PluginRegistryEntry {
  readonly structure?: {
    readonly themes?: PluginRegistryEntry.Themes;
    readonly representations?: readonly StructureRepresentationProvider<any>[];
    readonly presets?: {
      readonly hierarchy?: readonly TrajectoryHierarchyPresetProvider<any, any>[];
      readonly representation?: readonly StructureRepresentationPresetProvider<any, any>[];
    };
    readonly selectionQueries?: readonly StructureSelectionQuery[];
  };
  readonly volume?: {
    readonly themes?: PluginRegistryEntry.Themes;
    readonly representations?: readonly VolumeRepresentationProvider<any>[];
  };
  readonly particles?: {
    readonly themes?: PluginRegistryEntry.Themes;
    readonly representations?: readonly ParticleRepresentationProvider<any>[];
  };
  readonly formats?: readonly DataFormatProvider[];
  readonly lociLabels?: readonly LociLabelProvider[];
  readonly markdownExtensions?: readonly MarkdownExtension[];
  readonly dragAndDrop?: readonly PluginDragAndDropEntry[];
  readonly actions?: readonly (StateAction | StateTransformer | PluginSpec.Action)[];
  readonly animations?: readonly PluginStateAnimation[];
}

export namespace PluginRegistryEntry {
  export interface Themes {
    readonly color?: readonly ColorTheme.Provider<any, any>[];
    readonly size?: readonly SizeTheme.Provider<any, any>[];
  }
}

namespace PluginSpec {
  export interface Action {
    action: StateAction | StateTransformer;
    /* constructible react component with <action.customControl /> */
    customControl?: any;
    autoUpdate?: boolean;
  }

  export function Action(
    action: StateAction | StateTransformer,
    params?: {
      customControl?: any /* constructible react component with <action.customControl /> */;
      autoUpdate?: boolean;
    },
  ): Action {
    return { action, customControl: params && params.customControl, autoUpdate: params && params.autoUpdate };
  }

  export interface Behavior {
    transformer: StateTransformer;
    defaultParams?: any;
  }

  export function Behavior<T extends StateTransformer>(
    transformer: T,
    defaultParams: Partial<StateTransformer.Params<T>> = {},
  ): Behavior {
    return { transformer, defaultParams };
  }
}
