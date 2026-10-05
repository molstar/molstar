/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StateObject, StateObjectRef } from '@molstar/core/state';
import type { PluginContext } from '@molstar/plugin/context';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';

export interface PresetProvider<O extends StateObject = StateObject, P = any, S = {}> {
  id: string;
  display: { name: string; group?: string; description?: string };
  isApplicable?(a: O, plugin: PluginContext): boolean;
  params?(a: O | undefined, plugin: PluginContext): PD.For<P>;
  apply(a: StateObjectRef<O>, params: P, plugin: PluginContext): Promise<S> | S;
}
