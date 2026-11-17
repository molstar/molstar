/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultAssetManager } from '@molstar/core/util/assets';
import {
  type AssetPolicy,
  type AssetsPort,
  BrowserAssetPolicy,
  type SeqPorts,
  type TextMetricsPort,
} from '@molstar/seqx-core/context/index';
import { bundledWidth } from '@molstar/seqx-core/render';
import { browserTextMetrics } from '@molstar/seqx-viewer/canvas/index';

// The viewer's default ports (_spec/README.md "Ports"): each a factory with ownership, called once per viewer. Assets: a
// new `DefaultAssetManager` under the browser policy (https and relative locations under the page, no credentials,
// 256 MB an asset); text metrics: a canvas in the browser, the core's bundled table elsewhere. Both are disposed with
// the viewer; a host passes its own factories with `disposeWithViewer: false`.

/** The assets port's default adapter: the application's own manager is `{ manager, policy, disposeWithViewer: false }`. */
export function browserAssetsPort(policy: AssetPolicy = {}): () => AssetsPort {
  return () => ({
    manager: new DefaultAssetManager(),
    policy: {
      ...BrowserAssetPolicy,
      ...(typeof location !== 'undefined' ? { base: location.href } : {}),
      ...policy,
    },
  });
}

/** The text metrics port's default adapter: a canvas where there is a DOM, else the core's bundled table. */
export function defaultTextMetricsPort(): (() => TextMetricsPort) | undefined {
  return typeof document === 'undefined' ? undefined : () => browserTextMetrics() ?? { measure: bundledWidth };
}

/** The ports of the default spec. */
export function defaultPorts(policy: AssetPolicy = {}): SeqPorts {
  const textMetrics = defaultTextMetricsPort();
  return { assets: browserAssetsPort(policy), ...(textMetrics ? { textMetrics } : {}) };
}
