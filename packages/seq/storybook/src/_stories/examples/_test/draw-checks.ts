/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import type { FormatRuntime, OpenedSource, RegionRequest } from '@molstar/seqx-core/data';
import { type Frame, type RangesLay, recordDrawLists } from '@molstar/seqx-core/render';
import type { ExampleView } from '@molstar/seqx-examples';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { uiHistory } from '@molstar/seqx-viewer/ui-history';
import type { SeqViewerSpec } from '@molstar/seqx-viewer/viewer-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { readFileSync } from 'node:fs';
import { vi } from 'vitest';
import { withViewRegistry } from '../view-registry.js';

// The checks the example views' drawing tests share: every drawn mark of a ranges layer is hit (and outlined when
// hovered), and every item of the data in view has its marks.

export const HoverOwner = 'test';

export interface OpenOptions {
  /** The gene file's text: what every location outside `fixtures/` reads (the views that read the gene by location). */
  readonly gene?: string;
  /** Every location read from the examples' folder (the views of recorded datasets); `fixtures/…` always is. */
  readonly fixtures?: true;
  /** Each answer's delay, or `staggered` (the story's own: later regions later). */
  readonly delayMs?: number | 'staggered';
  /** How long the fake clock runs after the load. */
  readonly wait?: number;
  /** Keep the real `Date` and `performance` (a bound on time reads them): only the timers are faked. */
  readonly realClock?: boolean;
  /**
   * Await every runtime's open (a gzip file inflates on the platform's thread pool, by no clock) between turns of the
   * fake clock: what a view of recorded genome files needs to settle. `opens`: the opens alone, its region reads on the
   * fake clock (an answer's delay is a timer: awaited, it would wait for the clock the loop runs).
   */
  readonly awaitOpens?: true | 'opens';
  /** The host's history, for edits that are undone (the spec takes part in it). */
  readonly history?: boolean;
  /** The viewer's config items, as the story passes them (`exampleStory`'s `config`). */
  readonly config?: SeqViewerSpec['config'];
  /** Behaviors besides the default spec's, as the story lists them (`exampleStory`'s `behaviors`: `links.hover`). */
  readonly behaviors?: SeqViewerSpec['behaviors'];
  /** The viewer's extensions, as the story composes them (`exampleStory`'s `extensions`). */
  readonly extensions?: SeqViewerSpec['extensions'];
}

/** A viewer of an example view on a fake clock, composed as the stories compose it: loaded, the clock run, drawn. */
export async function openExample(view: ExampleView, options: OpenOptions): Promise<SeqViewer> {
  const { gene = '', delayMs = 300, wait = 6000, realClock = false, history = false } = options;
  vi.useFakeTimers({
    toFake: realClock ? ['setTimeout', 'clearTimeout'] : ['setTimeout', 'clearTimeout', 'Date', 'performance'],
  });
  const text = view.id === 'gene-by-region' || view.id === 'gene-region-ideogram' ? regionSeqText() : gene;
  const base = DefaultSeqViewerSpec();
  const opens: Promise<unknown>[] = [];
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    ...(history ? { history: uiHistory() } : {}),
    ...(options.config ? { config: [...(base.config ?? []), ...options.config] } : {}),
    ...(options.behaviors ? { behaviors: [...base.behaviors, ...options.behaviors] } : {}),
    ...(options.extensions ? { extensions: [...(base.extensions ?? []), ...options.extensions] } : {}),
    registry: withViewRegistry(
      base.registry ?? [],
      watchOpens(view.registry?.(delayMs === 'staggered' ? {} : { delayMs }) ?? [], opens, options.awaitOpens),
    ),
    ports: {
      assets: () => {
        // every dataset’s files (`fixtures/…`) are read as they are committed, at once (a read of the disk is not on the
        // fake clock: the layout waits for the letters); every other location answers the gene's text, or with
        // `fixtures` is read from the examples' folder too
        const read = (uri: string) => readFileSync(`${PackageFolder}${uri}`);
        const fixture = (uri: string) => options.fixtures === true || uri.startsWith('fixtures/');
        return {
          reader: {
            readText: async (uri) => (fixture(uri) ? read(uri).toString('utf8') : text),
            readBytes: async (uri) => (fixture(uri) ? new Uint8Array(read(uri)) : new TextEncoder().encode(text)),
          },
        };
      },
    },
  });
  v.resize(720);
  void v.load(view.build());
  if (options.awaitOpens) {
    await vi.advanceTimersByTimeAsync(200);
    for (let seen = -1; seen !== opens.length; ) {
      seen = opens.length;
      await Promise.allSettled(opens);
      await vi.advanceTimersByTimeAsync(200);
    }
  }
  // a lazy part the spec's model waits for (the square's chunk, a bridge's): its import is real I/O, never the fake
  // clock's, so it lands first (each load once)
  await Promise.resolve();
  for (
    let seen = new Set<unknown>(), m = v.spec.model.value;
    m.loading && !seen.has(m.loading);
    m = v.spec.model.value
  ) {
    seen.add(m.loading);
    await m.loading.catch(() => undefined);
  }
  await vi.advanceTimersByTimeAsync(wait);
  v.draw();
  return v;
}

/** The registry's runtimes with their opens kept (`awaitOpens`), else as they are. */
function watchOpens(
  entries: readonly SeqRegistryEntry[],
  opens: Promise<unknown>[],
  on: OpenOptions['awaitOpens'],
): SeqRegistryEntry[] {
  if (!on) return [...entries];
  return entries.map((e) =>
    e.kind === 'format-runtime'
      ? {
          ...e,
          runtime: {
            ...e.runtime,
            open: (source, ctx) => {
              // each region read too (not with `opens`): a settled view has none in flight
              const p = e.runtime.open(source, ctx).then((o) =>
                on === 'opens'
                  ? o
                  : Object.assign(Object.create(o) as OpenedSource, {
                      read: (request: RegionRequest, signal: AbortSignal) => {
                        const r = o.read(request, signal);
                        opens.push(r);
                        return r;
                      },
                    }),
              );
              opens.push(p);
              return p;
            },
          } satisfies FormatRuntime,
        }
      : e,
  );
}

/** The hover pass's rects. */
const hoverRects = (frame: Frame) =>
  recordDrawLists(frame.lists, ['hover'])
    .split('\n')
    .filter((l) => l.startsWith('rect ')).length;

/** Each drawn mark of every ranges layer hit at its centre: what is wrong, and how many marks were checked. */
export function checkMarks(v: SeqViewer) {
  const wrong: string[] = [];
  let marks = 0;
  let merged = 0;
  const layers = v.frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers));
  for (const l of layers) {
    if (l.type.name !== 'ranges') continue;
    const lay = l.lay as RangesLay;
    for (let j = 0; j < lay.m; j++) {
      const x = (lay.mx0[j] + lay.mx1[j]) / 2;
      const y = lay.laneY + (lay.mlane[j] + 0.5) * lay.laneHeight;
      // a lane scrolled half out of its band (or past a track not grown to its lanes yet) has its centre out of reach
      if (y < lay.box.y || y >= lay.box.y + lay.box.height) continue;
      marks++;
      // the items the mark draws: its own, or those of the run's marks
      const drawn = new Set<string | number>();
      if (lay.mitem[j] >= 0) drawn.add(lay.keys[lay.mitem[j]]);
      else {
        merged++;
        for (let k = 0; k < (lay.hn ?? 0); k++)
          if (lay.hlane[k] === lay.mlane[j] && lay.hx1[k] >= lay.mx0[j] && lay.hx0[k] <= lay.mx1[j])
            drawn.add(lay.keys[lay.hitem[k]]);
      }
      const hit = v.hit(x, y);
      const where = `${l.key} mark ${j} (lane ${lay.mlane[j]}, x ${lay.mx0[j]}–${lay.mx1[j]})`;
      if (hit?.layer !== l.key || hit.element === undefined || !drawn.has(hit.element.key)) {
        wrong.push(`${where}: ${hit?.layer} ${JSON.stringify(hit?.element)}`);
        continue;
      }
      // hovered, the item is outlined
      v.hover.set(hit, HoverOwner);
      const f = v.draw() ?? v.frame!;
      if (!hoverRects(f)) wrong.push(`${where}: hovered ${hit.element.key}, nothing outlined`);
      v.hover.set(undefined, HoverOwner);
    }
  }
  return { wrong, marks, merged };
}
