/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqViewerLayout, SeqViewerProvider } from '@molstar/seq-viewer/react';
import { renderSvg } from '@molstar/seq-viewer/svg';
import type { Meta, StoryObj } from './support/seqstar-story.js';
import { useCallback, useEffect, useMemo, useState } from 'react';
import { pDemoRuns, plasmid } from './support/seqstar-docs-views.js';
import f from './support/seqstar-format.module.css';
import { PrismCanvas, PrismPage } from './support/seqstar-prism.js';
import { useCaseContext, useStoryUi, ViewerCard } from './support/seqstar-viewer.js';
import s from './support/seqstar-viewer.module.css';
import { TextAction } from './ui/index.js';

/**
 * Cyclic coordinate systems: a sequence may be `cyclic`, its
 * positions taken modulo its length. A range on it may wrap the origin (start > stop), held as its two runs; a block
 * shows a piece across the origin (`range`) or the whole from another `origin`. Drawn linearly, and as a ring
 * with its ruler, marks, hover, keyboard cursor and an SVG export.
 */
const meta = { title: 'Seq*/Views and layouts/Axes and strands', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj;

export const Linear: Story = {
  name: 'A plasmid, linear',
  render: () => {
    const doc = useMemo(plasmid, []);
    return (
      <PrismPage
        title="A cyclic plasmid in linear views: a feature across the origin"
        caption={
          <>
            pDemo is 240 nucleotides and cyclic; its replication origin is written [200, 30): start after stop, across
            the origin. From the origin, it is drawn as its two runs, at the right end and the left, one item (hover
            either: both light) with no line between them. Rotated to origin 200, the ruler runs 201 … 240, 1 … 200 and
            the replication origin is one run at the start. A block's range may cross the origin too: [180, 60] shows
            181 … 240, 1 … 60. The last section is rotated to 150 and opens on the replication origin, its shortest arc.
            Lanes pack by arc: from the origin, the replication origin shares its lane with lacZα, which lies outside
            its arc. The same document as a ring: "A plasmid, ring".
          </>
        }
      >
        <ViewerCard
          doc={doc}
          label="Four linear sections over one cyclic sequence"
          options={{ parts: ['SectionHeader'] }}
        />
      </PrismPage>
    );
  },
};

/**
 * The ring, its marks set as a host would once the view is in (a selection and a highlight across the origin), and
 * beside it the same context exported by `renderSvg` with its overlay: again on request, with the hover as it stands.
 */
function RingAndExport() {
  const doc = useMemo(plasmid, []);
  const m = useCaseContext(doc, { view: 'ring' });
  const ui = useStoryUi(m);
  const [svg, setSvg] = useState<string | null>(null);
  const exported = useCallback(() => {
    if (m && !m.ctx.disposed)
      void renderSvg(m.ctx, { width: 480, overlay: true, hover: true, background: true }).then(setSvg);
  }, [m]);
  useEffect(() => {
    if (!m) return;
    let done = false;
    const sub = m.ctx.graph.assembled.subscribe((a) => {
      if (done || !a?.view || m.ctx.views.current.value !== 'ring') return;
      done = true;
      // The selection 226 … 12 and the highlight 236 … 6, both across the origin; then the export of the figure.
      m.ctx.selection.set(pDemoRuns(0, 225, 240, 0, 0, 12), 'api');
      m.ctx.highlight.set('story', pDemoRuns(0, 235, 240, 0, 0, 6));
      exported();
    });
    return () => sub.unsubscribe();
  }, [m, exported]);
  return (
    <div className={s.pair}>
      <section className={f.card}>
        <h2 className={f.heading}>ring: one circular section over the cyclic sequence</h2>
        {m ? (
          <SeqViewerProvider context={m.ctx} ui={ui}>
            <PrismCanvas>
              <SeqViewerLayout />
            </PrismCanvas>
          </SeqViewerProvider>
        ) : (
          <div className={f.status}>Loading…</div>
        )}
      </section>
      <section className={f.card}>
        <div className={f.row}>
          <h2 className={f.heading}>renderSvg(ctx, {'{ width: 480, overlay: true, hover: true }'})</h2>
          <TextAction tone="quiet" onClick={exported}>
            Export again
          </TextAction>
        </div>
        {svg ? (
          <div className={s.figure} dangerouslySetInnerHTML={{ __html: svg }} />
        ) : (
          <div className={f.status}>Drawing…</div>
        )}
      </section>
    </div>
  );
}

export const Ring: Story = {
  name: 'A plasmid, ring',
  render: () => (
    <PrismPage
      title="A cyclic plasmid as a ring: the circular geometry"
      caption={
        <>
          The same tracks with the section's style "circular": rings from the outside in (the ruler, the letters as
          cells, the features, the cut sites), the axis from the top and clockwise. The ruler's numbers sit outside the
          ring, with unlabeled ticks every 2 between them. The replication origin, [200, 30), crosses the top as one arc
          on the features' ring, sharing it with lacZα. A selection, 226 … 12, and a highlight, 236 … 6, both cross the
          origin: each is one sector over every ring, from the ruler's ticks to the hub, the selection with its edges.
          Hover a position: a radial line; a feature: its sector, and its outline on its ring. Tab into the ring and use
          the arrows: the keyboard cursor outlines its cell's sector and goes round past 240 onto 1 (Home: 1, End: 240;
          up and down change rings). Beside it, the same context exported as SVG with its overlay (the marks, and the
          hover as it stands when exported). Labels around the ring and chords are not drawn yet.
        </>
      }
    >
      <RingAndExport />
    </PrismPage>
  ),
};
