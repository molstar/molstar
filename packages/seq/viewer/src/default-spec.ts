/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqRuntime } from '@molstar/seqx-core/data';
import {
  AddressingProvider,
  CursorProvider,
  EditingProvider,
  HoverProvider,
  LegendProvider,
  NavigateProvider,
  NavigatorProvider,
  SelectProvider,
  RowsProvider,
  WindowProvider,
} from '@molstar/seqx-viewer/behaviors';
import { defaultPorts } from '@molstar/seqx-viewer/ports';
import { frameAnimation } from './animation.js';
import { legendsMount } from './mount-legend.js';
import { navigatorsMount } from './mount-navigator.js';
import type { SeqViewerSpec } from './viewer-spec.js';

// The default spec (_spec/README.md "Composition"): the built-in entries, behaviors, ports and UI, as data. An
// application spreads it and changes what it needs, or writes a smaller spec; what is not listed is not in its bundle.
// Each call makes a spec of its own (its ports are factories called once per viewer). The history is the UI's: a
// mounted viewer gets a `@molstar/ui/history` stack of its own (`ui.history`), a host gives its stack (`history`,
// `hostHistory`), and a headless spread (`ui: undefined`) has none. The clock is the page's frames (`frameAnimation`): a
// host with a clock of its own, or a headless export, gives its adapter (`animation`).

export function DefaultSeqViewerSpec(): SeqViewerSpec {
  return {
    registry: [{ kind: 'format-runtime', runtime: SeqRuntime }],
    behaviors: [
      { provider: HoverProvider },
      { provider: NavigateProvider },
      { provider: SelectProvider },
      { provider: CursorProvider },
      { provider: AddressingProvider },
      { provider: EditingProvider },
      { provider: LegendProvider },
      { provider: WindowProvider },
      { provider: NavigatorProvider },
      { provider: RowsProvider },
    ],
    bindings: {},
    config: [],
    ports: defaultPorts(),
    animation: frameAnimation(),
    ui: { canvas: true, history: true, legends: legendsMount, navigators: navigatorsMount },
  };
}
