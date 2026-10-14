/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Loaders } from '@molstar/seq-core/format';
import { SeqCommands } from '@molstar/seq-viewer';
import { LETTERS_PER_PX } from '@molstar/seq-viewer/draw';
import { useLiveVersion, useSeqViewer } from '@molstar/seq-viewer/react';
import type { Meta, StoryObj } from '@storybook/react-vite';
import { useMemo } from 'react';
import { lazyStoryDoc, lazyStoryLoaders } from '../../../viewer/_test/lazy-story-fixtures.js';
import {
  STREAMED_AT as AT,
  CODING_PROTEIN,
  codingLoaders,
  codingStoryDoc,
  late,
  STREAMED_SECTION as SECTION,
  streamedDoc,
  streamedKnobs,
  streamedLoaders,
  transcriptLoaders,
  transcriptStoryDoc,
} from './support/seqstar-docs-streamed.js';
import f from './support/seqstar-format.module.css';
import { PrismPage } from './support/seqstar-prism.js';
import { ViewerCard } from './support/seqstar-viewer.js';
import { TextAction } from './ui/index.js';

/**
 * Streamed data: a chromosome whose letters a stream serves by region,
 * as a 2bit file would, from a committed fixture provider (offline: synthetic letters, N runs, soft-masked stretches).
 * Nothing is asked until a device pixel holds at most `LETTERS_PER_PX` positions; then the regions in view arrive (4,096
 * letters each, the view's first, then a view width each side) and are drawn as cells, and as letters once they fit.
 * Genes (exons of one id joined, a lane each kept across regions), variants (by alternative) and a values track
 * from the same provider; zoomed out, genes and variants as density from `count` bins. A region asked and not yet
 * arrived shows a quiet placeholder band (the latency knob makes it last); past the letters gate, a hatched band and a
 * hint. Progressive loading, on the overlay: skeletons that shimmer, error bands (the failure knob) with Retry, a
 * fade as regions arrive; genes colored by category (the manifest's vocabulary); and a lazy alignment loaded by view.
 */
const meta = { title: 'Seq*/Streaming/Streamed data', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj<{ latency: number; failure: number }>;

const ZOOMS: [string, number][] = [
  ['2 Mb', 2_000_000],
  ['100 kb', 100_000],
  ['10 kb', 10_000],
  ['1 kb', 1000],
  ['100 positions', 100],
];

/** The zoom buttons, Retry, and what the stream manager has done. */
function Panel({ zooms = true }: { zooms?: boolean }) {
  const { context: ctx } = useSeqViewer();
  useLiveVersion(ctx.streams.stats);
  const st = ctx.streams.stats.current;
  const retry = () => void ctx.commands.run({ command: SeqCommands.Retry.id, params: {} }, 'api');
  const go = (span: number) => {
    const vp = ctx.viewports.get(SECTION);
    if (span === 0) vp?.fit();
    else vp?.zoomTo(0, AT, AT + span);
  };
  return (
    <div className={f.stack}>
      <div className={f.row}>
        {zooms && <TextAction onClick={() => go(0)}>Whole chromosome</TextAction>}
        {zooms &&
          ZOOMS.map(([label, span]) => (
            <TextAction key={label} onClick={() => go(span)}>
              {label}
            </TextAction>
          ))}
        <TextAction onClick={retry}>Retry failed regions (seqv:retry)</TextAction>
      </div>
      <div className={f.code}>
        {st.requests} regions asked, {st.inFlight} in flight, {st.cancelled} cancelled, {st.errors} failed;{' '}
        {(st.bytes / 1024).toFixed(0)} KiB held
      </div>
    </div>
  );
}

function Card({ latency, failure }: { latency: number; failure: number }) {
  streamedKnobs.latency = latency;
  streamedKnobs.failure = failure;
  const d = useMemo(streamedDoc, []);
  const options = useMemo(() => ({ loaders: streamedLoaders(), parts: ['TrackActions'] as const }), []);
  return (
    <ViewerCard doc={d} options={options} label="A 20 Mb chromosome, everything streamed by region">
      <Panel />
    </ViewerCard>
  );
}

const argTypes = {
  latency: { control: { type: 'range', min: 0, max: 3000, step: 50 }, description: "ms before each region's answer" },
  failure: {
    control: { type: 'range', min: 0, max: 1, step: 0.05 },
    description: 'the share of region requests that fail (Retry asks them again)',
  },
} as const;

export const Letters: Story = {
  name: 'Letters, genes, variants and a signal',
  args: { latency: 250, failure: 0 },
  argTypes,
  render: (args) => (
    <PrismPage
      title="Streamed data: a chromosome's letters, genes, variants and signal, by region"
      caption={
        <>
          The letters come from a stream, as from a 2bit file: the manifest gives the chromosome's length, and nothing
          else is asked until you zoom in to {LETTERS_PER_PX} positions a device pixel or fewer. Then the regions in
          view arrive (4,096 letters each, the view's first, then a view width on each side), drawn as cells colored by
          nucleotide, and as letters once a position is wide enough to read one. Pan: the margin's regions are already
          there. Every 100 kb starts with a run of N; soft-masked stretches arrive in lower case and are drawn as upper
          case. Zoomed out, the Sequence track is hatched and says to zoom in.
        </>
      }
    >
      <div className={f.kicker}>
        Genes and variants come from items streams: zoomed out (2 Mb and wider), as density from each region's counts
        per bin; closer, as items, 65,536 positions a region. A gene's exons are joined, and it keeps its lane and is
        drawn once though it crosses regions, colored by its category (protein coding, lncRNA, pseudogene: the
        manifest's vocabulary gives the palette its order); hover one for its label and fields. Variants are colored by
        their alternative, a second one at a position stacked above.
      </div>
      <div className={f.kicker}>
        Progressive loading: a region asked and not yet arrived is a skeleton, the pending shade with a slow shimmer
        (one sweep a second); the signal keeps its coarser level drawn instead. Raise the latency and pan or jump:
        regions fade in over 150 ms as they arrive, those arriving together together. Raise the failure share: failed
        regions are hatched red bands saying so, the hover label says "failed", and the track's "…" menu (or the button,
        `seqv:retry`) asks them again. With reduced motion the shade holds still.
      </div>
      <Card latency={args.latency} failure={args.failure} />
    </PrismPage>
  ),
};

// ---- A protein on its genome (ST7) ---------------------------------------------------------------------------

function CodingCard({ latency, failure }: { latency: number; failure: number }) {
  streamedKnobs.latency = latency;
  streamedKnobs.failure = failure;
  const d = useMemo(codingStoryDoc, []);
  const options = useMemo(() => ({ loaders: codingLoaders(), parts: ['TrackActions'] as const }), []);
  return (
    <ViewerCard doc={d} options={options} label="P1, 1,334 residues coded by three exons over 220 kb of chr2">
      <CodingZooms />
    </ViewerCard>
  );
}
function CodingZooms() {
  const { context: ctx } = useSeqViewer();
  useLiveVersion(ctx.streams.stats);
  const st = ctx.streams.stats.current;
  const go = (a: number, e: number) => ctx.viewports.get(CODING_PROTEIN)?.zoomTo(0, a, e);
  return (
    <div className={f.stack}>
      <div className={f.row}>
        <TextAction onClick={() => go(0, 1334)}>Whole protein</TextAction>
        <TextAction onClick={() => go(450, 550)}>Residues 450–550 (a codon split at 500)</TextAction>
        <TextAction onClick={() => go(880, 920)}>Residues 880–920</TextAction>
      </div>
      <div className={f.code}>
        {st.requests} regions asked, {st.inFlight} in flight; {(st.bytes / 1024).toFixed(0)} KiB held
      </div>
    </div>
  );
}

export const ThroughAPath: Story = {
  name: 'A protein on its genome',
  args: { latency: 250, failure: 0 },
  argTypes,
  render: (args) => (
    <PrismPage
      title="Streamed data through a path: a genome's tracks on a protein's axis"
      caption={
        <>
          P1 is a protein whose gene sits on chr2 of a streamed genome: three exons over 220 kb, joined to the protein
          by a codon mapping. Its section shows chr2's streamed genes, variants and phyloP on the protein's own axis,
          through that mapping: only the regions its exons reach are asked (three regions of genes, never the introns
          between them), a variant lands on the residue its position codes, a gene's exon on its coding part, and phyloP
          over a residue is read over its codon (both parts where an intron splits one, as at residues 500 and 900).
          Hover a variant or a gene: its range is lit on the protein. Hover any residue: its label shows the codon that
          codes it, read from the genome's streamed letters (asked when first hovered), with its positions on chr2 and
          its strand; a codon an intron splits shows both parts (residue 500's T|CA, 900's A|TG). The second section is
          chr2 at the gene, the same tracks shown directly.
        </>
      }
    >
      <CodingCard latency={args.latency} failure={args.failure} />
    </PrismPage>
  ),
};

// ---- A protein on its streamed transcript (ST7) ---------------------------------------------------------------

function TranscriptCard({ latency, failure }: { latency: number; failure: number }) {
  streamedKnobs.latency = latency;
  streamedKnobs.failure = failure;
  const d = useMemo(transcriptStoryDoc, []);
  const options = useMemo(() => ({ loaders: transcriptLoaders(), parts: ['TrackActions'] as const }), []);
  return <ViewerCard doc={d} options={options} label="T1, a streamed transcript on chr2, and P1 mapped onto it" />;
}

export const OnItsTranscript: Story = {
  name: 'A protein on its streamed transcript',
  args: { latency: 250, failure: 0 },
  argTypes,
  render: (args) => (
    <PrismPage
      title="Item coordinate systems: a protein's domains on the streamed transcript that codes it"
      caption={
        <>
          The genes stream serves transcript T1 with a map, its residues over its CDS on chr2 (a codon mapping's
          params). When a region holding it arrives, T1 becomes a coordinate system of its own, genes_src/T1, mapped
          onto chr2. The document maps protein P1 onto it, base for base, before it exists: P1's domains appear along
          T1's exons in the genome section as soon as T1 does, a domain across an intron in two parts. The second
          section is T1's own residues (its block names genes_src/T1): P1's domains, and the genome's variants and
          phyloP through T1's map. Hover a residue there: it is numbered in T1's own positions, and its codon is read
          from the genome's letters.
        </>
      }
    >
      <TranscriptCard latency={args.latency} failure={args.failure} />
    </PrismPage>
  ),
};

// ---- A lazy alignment by view ---------------------------------------------------------------------------------

/** The story's document and transport (shared with viewer/_test/lazy-story.test.ts, which opens them headlessly). */
const lazyLoaders: Loaders = lazyStoryLoaders((fetch) => late(fetch));
function LazyCard({ latency, failure }: { latency: number; failure: number }) {
  streamedKnobs.latency = latency;
  streamedKnobs.failure = failure;
  const d = useMemo(lazyStoryDoc, []);
  const options = useMemo(() => ({ loaders: lazyLoaders, maxHeight: 460, parts: ['TrackActions'] as const }), []);
  return (
    <ViewerCard doc={d} options={options} label="20,000 rows × 1,200 columns, tiles loaded by view">
      <Panel zooms={false} />
    </ViewerCard>
  );
}

export const LazyAlignment: Story = {
  name: 'A lazy alignment by view',
  args: { latency: 600, failure: 0 },
  argTypes,
  render: (args) => (
    <PrismPage
      title="Streamed data: a lazy alignment, tile by tile"
      caption={
        <>
          A family alignment too large to hold, served in 256 × 256 tiles at the level the view needs (the reference
          protocol, `seqstar-tiles/1`, over a sampled transport: letters from a hash of the cell). The coarsest level
          comes first and stays drawn under a finer tile until it arrives; a tile with nothing in its place is a
          skeleton that shimmers, and arriving tiles fade in. Raise the latency and scroll or zoom to watch; raise the
          failure share for hatched red error bands, then Retry from the Rows track's "…" menu or the button.
        </>
      }
    >
      <LazyCard latency={args.latency} failure={args.failure} />
    </PrismPage>
  ),
};
