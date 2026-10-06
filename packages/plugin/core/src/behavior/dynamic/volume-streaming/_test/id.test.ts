/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject } from '@molstar/plugin/state/objects';
import { VolumeStreaming } from '../behavior.js';
import { VolumeStreamingObjectType, isVolumeStreamingObject } from '../id.js';

describe('volume streaming state object id', () => {
  it('is the type of the VolumeStreaming state object', () => {
    expect(VolumeStreaming.type).toBe(VolumeStreamingObjectType);
    expect(VolumeStreamingObjectType).toEqual({ name: 'Volume Streaming', typeClass: 'Behavior' });
  });

  it('recognizes VolumeStreaming objects without the behavior', () => {
    const obj = new VolumeStreaming({} as VolumeStreaming.Behavior);
    expect(isVolumeStreamingObject(obj)).toBe(true);
    expect(isVolumeStreamingObject(obj)).toBe(VolumeStreaming.is(obj));
  });

  it('rejects other objects', () => {
    const other = new PluginStateObject.Group({});
    expect(isVolumeStreamingObject(other)).toBe(false);
    expect(isVolumeStreamingObject(undefined)).toBe(false);
  });
});
