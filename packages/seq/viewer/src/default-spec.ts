/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqRuntime } from '@molstar/seqx-core/data';
import { HoverProvider, NavigateProvider } from '@molstar/seqx-viewer/behaviors';
import { defaultPorts } from '@molstar/seqx-viewer/ports';
import type { SeqViewerSpec } from './viewer-spec.js';

// The default spec (_spec/README.md "Composition"): the built-in entries, behaviors, ports and UI, as data. An
// application spreads it and changes what it needs, or writes a smaller spec; what is not listed is not in its bundle.
// Each call makes a spec of its own (its ports are factories called once per viewer).

export function DefaultSeqViewerSpec(): SeqViewerSpec {
  return {
    registry: [{ kind: 'format-runtime', runtime: SeqRuntime }],
    behaviors: [{ provider: HoverProvider }, { provider: NavigateProvider }],
    bindings: {},
    config: [],
    ports: defaultPorts(),
    ui: { canvas: true },
  };
}
