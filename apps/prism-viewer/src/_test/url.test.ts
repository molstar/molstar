/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { resolveOptions } from '@molstar/prism-viewer/options';
import { modesFromUrl, optionsFromUrl, unknownUrlParameters } from '../url.js';

describe('unknownUrlParameters', () => {
  const load = ['pdb', 'url', 'mvs-url'];

  it('names the parameters neither list knows, once each', () => {
    expect(unknownUrlParameters('?pdb=1cbs&hide-controls=1&foo=1&bar&foo=2&PIXEL-SCALE=2', load)).toEqual([
      'foo',
      'bar',
    ]);
  });

  it('knows the options, the dropped ones and the load parameters', () => {
    expect(unknownUrlParameters('?pick-padding=1&map-provider=rcsb&mvs-url=x&debug-mode=1', load)).toEqual([]);
    expect(unknownUrlParameters('', load)).toEqual([]);
  });

  it('knows the structure load options', () => {
    expect(unknownUrlParameters('?pdb=1cbs&assembly=1&all-models=1', load)).toEqual([]);
  });
});

describe('optionsFromUrl', () => {
  it("maps v5's map-provider to the volume streaming server; an explicit server wins", () => {
    expect(optionsFromUrl('?map-provider=rcsb').data?.volumeStreamingServer).toBe('https://maps.rcsb.org');
    expect(optionsFromUrl('?map-provider=pdbe').data?.volumeStreamingServer).toBe(
      'https://www.ebi.ac.uk/pdbe/densities',
    );
    expect(optionsFromUrl('?map-provider=nope').data?.volumeStreamingServer).toBeUndefined();
    expect(optionsFromUrl('?map-provider=rcsb&volume-streaming-server=https://x.org').data?.volumeStreamingServer).toBe(
      'https://x.org',
    );
  });

  it('sets the fullscreen button only when show-toggle-fullscreen is given, keeping the other buttons', () => {
    expect(optionsFromUrl('').ui?.buttons).toBeUndefined();
    expect(optionsFromUrl('?show-toggle-fullscreen=0').ui?.buttons).toEqual({ fullscreen: false });
    const options = resolveOptions(optionsFromUrl('?show-toggle-fullscreen=1'));
    expect(options.ui.buttons.fullscreen).toBe(true);
    expect(options.ui.buttons.reset).toBe(true);
    expect(resolveOptions(optionsFromUrl('?show-toggle-fullscreen=0')).ui.buttons.fullscreen).toBe(false);
  });

  it('reads the debug and timing modes', () => {
    expect(modesFromUrl('?debug-mode=1&timing-mode=1')).toEqual({ debug: true, timing: true });
    expect(modesFromUrl('?debug-mode=0')).toEqual({ debug: false, timing: false });
  });
});
