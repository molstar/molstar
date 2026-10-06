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

export { PluginSpec };

interface PluginSpec {
  actions?: PluginSpec.Action[];
  behaviors: PluginSpec.Behavior[];
  animations?: PluginStateAnimation[];
  customFormats?: [string, DataFormatProvider][];
  canvas3d?: PartialCanvas3DProps;
  layout?: {
    initial?: Partial<PluginLayoutStateProps>;
  };
  config?: [PluginConfigItem, unknown][];
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
