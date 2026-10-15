/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';
import { MolstarViewerExtensions } from './extensions.js';

/**
 * The viewer with Mol*'s extensions (app.md §2): `DefaultPrismViewerSpec()` (the built-ins) and the app's extension
 * list, without the assistant (a host that wants it uses the app, or composes the packages). The name differs from
 * the default spec's on purpose.
 */
export function FullPrismViewerSpec(): PrismViewerSpec {
  const base = DefaultPrismViewerSpec();
  return { ...base, extensions: [...(base.extensions ?? []), ...MolstarViewerExtensions] };
}
