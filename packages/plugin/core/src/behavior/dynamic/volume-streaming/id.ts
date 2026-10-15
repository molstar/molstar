/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StateObject } from '@molstar/core/state';
import type { PluginStateObject } from '@molstar/plugin/state/objects';

/**
 * State object type of the volume streaming behavior. Shared between the `VolumeStreaming` state object (in
 * `behavior.ts`) and base modules that only need to recognize it, so they do not load the behavior.
 */
export const VolumeStreamingObjectType: PluginStateObject.TypeInfo = {
  name: 'Volume Streaming',
  typeClass: 'Behavior',
};

/** Checks whether a state object is the volume streaming behavior object. */
export function isVolumeStreamingObject(obj?: StateObject): boolean {
  return !!obj && obj.type === VolumeStreamingObjectType;
}
