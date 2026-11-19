/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Size entry 4, a minimal viewer: composed with tracks, letters and ranges only (the seq runtime, hover and navigation;
// no navigators, stories, editing, history, records runtime, rows, grids). Bundled by `sizes.mjs`, never run; the
// script checks its bundle holds none of what it leaves out (`Unlisted` there).

import { SeqRuntime } from '@molstar/seqx-core/data';
import { LettersType, RangesType, TracksLayout } from '@molstar/seqx-core/render';
import { createSeqViewer } from '@molstar/seqx-viewer';
import { frameAnimation } from '@molstar/seqx-viewer/animation';
import { HoverProvider, NavigateProvider } from '@molstar/seqx-viewer/behaviors';
import { defaultPorts } from '@molstar/seqx-viewer/ports';

export const create = (target: HTMLElement) =>
  createSeqViewer(
    {
      registry: [
        { kind: 'format-runtime', runtime: SeqRuntime },
        { kind: 'layer-type', type: LettersType },
        { kind: 'layer-type', type: RangesType },
        { kind: 'layout', layout: TracksLayout },
      ],
      behaviors: [{ provider: HoverProvider }, { provider: NavigateProvider }],
      ports: defaultPorts(),
      animation: frameAnimation(),
      ui: { canvas: true },
    },
    target,
  );
