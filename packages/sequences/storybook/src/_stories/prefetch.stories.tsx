/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqViewerConfig } from '@molstar/seq-viewer';
import { useSeqViewer } from '@molstar/seq-viewer/react';
import type { Meta, StoryObj } from './support/seqstar-story.js';
import { useEffect, useMemo, useSyncExternalStore } from 'react';
import f from './support/seqstar-format.module.css';
import {
  AskedLog,
  PREFETCH_SECTION,
  PREFETCH_SNAPSHOTS,
  prefetchDoc,
  prefetchLoaders,
} from './support/seqstar-prefetch.js';
import { PrismPage } from './support/seqstar-prism.js';
import { ViewerCard } from './support/seqstar-viewer.js';
import { ScrollArea, TextAction } from './ui/index.js';

/**
 * Prefetch: the viewer asks for regions just ahead of where the view is going, so data is usually there when it
 * arrives. Beside the view is the log of every request the providers get, in order, the prefetches marked.
 */
const meta = { title: 'Seq*/Streaming/Prefetch', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj<{ prefetch: boolean }>;

const mb = (n: number) => (n / 1e6).toFixed(2);
/** The snapshot the document opens on (its focus is the opening view). */
const [first] = PREFETCH_SNAPSHOTS;

/** The request log: newest first, a prefetch marked, a coarse region (a level) by its level. */
function Log({ log }: { log: AskedLog }) {
  const entries = useSyncExternalStore(log.subscribe, log.snapshot);
  const shown = entries.slice(-80).reverse();
  const ahead = entries.filter((e) => e.priority === 'prefetch').length;
  return (
    <div className={f.stack}>
      <div className={f.row}>
        <span className={f.code}>
          {entries.length} requests, {ahead} prefetches
        </span>
        <TextAction onClick={() => log.clear()}>Clear</TextAction>
      </div>
      <ScrollArea className={f.box} orientation="vertical">
        <table className={f.table}>
          <thead>
            <tr>
              <th>#</th>
              <th>asked</th>
              <th>source</th>
              <th>level</th>
              <th>Mb</th>
            </tr>
          </thead>
          <tbody>
            {shown.map((e) => (
              <tr key={e.n} style={e.priority === 'prefetch' ? { fontWeight: 600 } : undefined}>
                <td className={f.mono}>{e.n}</td>
                <td>{e.priority === 'prefetch' ? 'prefetch ›' : e.priority}</td>
                <td>{e.source}</td>
                <td className={f.mono}>{e.level < 0 ? 'data' : e.level}</td>
                <td className={f.mono}>
                  {mb(e.start)}–{mb(e.stop)}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </ScrollArea>
    </div>
  );
}

/** What moves the view: shifts to either side, zooms, and the snapshots (entered from the cache, if hinted). */
function Panel({ prefetch }: { prefetch: boolean }) {
  const { context: ctx } = useSeqViewer();
  const vp = () => ctx.viewports.get(PREFETCH_SECTION);
  useEffect(() => ctx.config.set(SeqViewerConfig.Stream.Prefetch, prefetch), [ctx, prefetch]);
  const pan = (px: number) => vp()?.panBy(px, 0, 'api');
  const zoom = (k: number) => {
    const v = ctx.layout.axis(PREFETCH_SECTION);
    const [a, b] = [v?.visLo[0] ?? 0, v?.visHi[0] ?? 0];
    const mid = (a + b) / 2;
    vp()?.zoomTo(0, Math.round(mid - ((b - a) * k) / 2), Math.round(mid + ((b - a) * k) / 2));
  };
  return (
    <div className={f.stack}>
      <div className={f.row}>
        <TextAction onClick={() => vp()?.zoomTo(0, first[1], first[2])}>
          Where it opens: 100 kb at {mb(first[1])} Mb
        </TextAction>
        <TextAction onClick={() => pan(-300)}>‹ Pan left</TextAction>
        <TextAction onClick={() => pan(300)}>Pan right ›</TextAction>
        <TextAction onClick={() => zoom(0.5)}>Zoom in</TextAction>
        <TextAction onClick={() => zoom(2)}>Zoom out</TextAction>
      </div>
      <div className={f.row}>
        {PREFETCH_SNAPSHOTS.map(([title]) => (
          <TextAction key={title} onClick={() => void ctx.snapshots.apply(title)}>
            Snapshot {title}
          </TextAction>
        ))}
        <TextAction onClick={() => void ctx.snapshots.apply(null)}>Base</TextAction>
      </div>
    </div>
  );
}

function Card({ prefetch }: { prefetch: boolean }) {
  const log = useMemo(() => new AskedLog(), []);
  const doc = useMemo(prefetchDoc, []);
  const options = useMemo(() => ({ loaders: prefetchLoaders(log) }), [log]);
  return (
    <div style={{ display: 'grid', gridTemplateColumns: 'minmax(0, 2fr) minmax(280px, 1fr)', gap: 16 }}>
      <ViewerCard doc={doc} options={options} label="A 20 Mb chromosome, three snapshots">
        <Panel prefetch={prefetch} />
      </ViewerCard>
      <section className={f.card}>
        <h2 className={f.heading}>Requests</h2>
        <Log log={log} />
      </section>
    </div>
  );
}

export const Ahead: Story = {
  name: 'Requests ahead of the view',
  args: { prefetch: true },
  argTypes: {
    prefetch: { control: 'boolean', description: "`Stream.Prefetch`: off, only the view's and the margin's are asked" },
  },
  render: (args) => (
    <PrismPage
      title="Prefetch: regions asked just ahead of where the view is going"
      caption={
        <>
          Pan the view and the log shows, after the view's own regions and the margin's, a view width more in the
          direction it went (a prefetch, in bold): it is kept while the view rests and moves to the other side when you
          pan back. Zoom in and the next finer level of what is in view is asked; zoom out, the next coarser. The
          snapshots A, B and C focus 100 kb each far apart: entering one asks for the regions of its neighbors, so the
          next one is drawn from what is held when you get there. A prefetch is never waited for (the view settles
          without it), and a view's region takes the slot of a prefetch when the in-flight bound is full. Turn the knob
          off to see only the view's and the margin's.
        </>
      }
    >
      <Card prefetch={args.prefetch} />
    </PrismPage>
  ),
};
