/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import type { ExampleView } from '@molstar/seqx-examples';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import type { SeqViewerSpec } from '@molstar/seqx-viewer/viewer-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { vi } from 'vitest';
import type { ChainTable } from '../../../support/fake-structure/chain-table.js';
import { chainTableOf } from '../compose.js';

// What the structure link's twins share: a headless viewer of an example view (its files read from the examples'
// fixtures, as the page reads them), the dataset's chain table from disk, and the fake clock the stories' timers run on.

export const readFixture = (file: string) => readFileSync(`${PackageFolder}${file}`, 'utf8');
export const tableOf = (file: string): ChainTable => chainTableOf(readFixture(file));

/** A settled headless viewer of the view (real timers: call `fakeClock` after). */
export async function openView(view: ExampleView, spec: Partial<SeqViewerSpec> = {}): Promise<SeqViewer> {
  const base = DefaultSeqViewerSpec();
  const read = (uri: string) => readFileSync(`${PackageFolder}${uri}`);
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    ...spec,
    registry: [...(base.registry ?? []), ...(spec.registry ?? [])],
    behaviors: [...base.behaviors, ...(spec.behaviors ?? [])],
    ports: {
      assets: () => ({
        reader: {
          readText: async (uri) => read(uri).toString('utf8'),
          readBytes: async (uri) => new Uint8Array(read(uri)),
        },
      }),
    },
  });
  v.resize(1440);
  await v.load(view.build());
  await v.settle();
  return v;
}

/** The stories' clock: timers faked, a frame a timer 16 ms on; returns the restore. */
export function fakeClock(): () => void {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const raf = DomHost.requestAnimationFrame;
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  return () => {
    DomHost.requestAnimationFrame = raf;
    vi.useRealTimers();
  };
}
