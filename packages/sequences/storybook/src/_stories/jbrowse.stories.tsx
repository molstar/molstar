/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Loaders, StreamProvider } from '@molstar/seq-core/format';
import {
  bigBedProvider,
  bigWigProvider,
  chrNames,
  fetchRanges,
  renamed,
  twoBitProvider,
  vcfProvider,
} from '@molstar/seq-core/stream/gmod';
import { SeqCommands, SeqViewerConfig } from '@molstar/seq-viewer';
import { useLiveVersion, useSeqViewer } from '@molstar/seq-viewer/react';
import type { Meta, StoryObj } from './support/seqstar-story.js';
import { useEffect, useMemo, useState } from 'react';
import { GRCH38 } from '../../../core/_test/format/stream-provider.js';
import { JBROWSE_SECTION, jbrowseDoc } from '../../../viewer/_test/jbrowse-story-fixtures.js';
import f from './support/seqstar-format.module.css';
import { IDEOGRAM_REGISTRY } from './support/seqstar-ideogram.js';
import { PrismPage } from './support/seqstar-prism.js';
import { ViewerCard } from './support/seqstar-viewer.js';
import { TextAction } from './ui/index.js';

/**
 * The genome browser demo: hg38 as JBrowse 2's own demo serves it, read **live**
 * from the public files through `@molstar/seq-core/stream/gmod` over HTTP ranges (the one story that needs the network): the
 * genome's letters from UCSC's 2bit, GENCODE's genes from its BigBed, GC content from a BigWig and ClinVar from its
 * tabix VCF (named `1`, `MT`: mapped to `chr1`, `chrM` by `renamed` and `chrNames`). Everything is asked by view, the
 * VCF's dense regions through the zoom gate (counts, then finer regions where a pixel holds 16 positions or fewer).
 */
const meta = { title: 'Seq*/Cases/Genome browser', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj<{ shimmerFps: number }>;

const UCSC = 'https://hgdownload.soe.ucsc.edu';
const FILES = {
  genome: `${UCSC}/goldenPath/hg38/bigZips/hg38.2bit`,
  gc: `${UCSC}/gbdb/hg38/bbi/gc5BaseBw/gc5Base.bw`,
  genes: `${UCSC}/gbdb/hg38/gencode/gencodeV50.bb`,
  clinvar: 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2024/clinvar_20240107.vcf.gz',
};
/** The primary assembly's chromosomes and their lengths (the VCF's index has none), chrM's from the 2bit. */
const LENGTHS = new Map<string, number>([...GRCH38.map(([k, n]) => [k, n] as [string, number]), ['chrM', 16_569]]);
const primary = (name: string) => (LENGTHS.has(name) ? name : null);
const CLASSES = ['coding', 'nonCoding', 'pseudo'];

/** The four providers, made once per page (each keeps its own small caches of the file it reads). */
let made: Record<string, StreamProvider> | null = null;
function providers(): Record<string, StreamProvider> {
  made ??= {
    genome: twoBitProvider(fetchRanges(FILES.genome), { only: (n) => LENGTHS.has(n) }),
    gc: renamed(bigWigProvider(fetchRanges(FILES.gc)), primary),
    genes: renamed(
      bigBedProvider(fetchRanges(FILES.genes), { label: 'geneName', category: 'transcriptClass', vocabulary: CLASSES }),
      primary,
    ),
    clinvar: renamed(
      vcfProvider(
        fetchRanges(FILES.clinvar),
        fetch(`${FILES.clinvar}.tbi`).then(async (r) => new Uint8Array(await r.arrayBuffer())),
        {
          info: ['CLNSIG', 'CLNDN', 'GENEINFO'],
          lengths: Object.fromEntries([...LENGTHS].map(([k, n]) => [k === 'chrM' ? 'MT' : k.slice(3), n])),
        },
      ),
      chrNames,
    ),
  };
  return made;
}
const loaders: Loaders = {
  load: () => Promise.reject(new Error('no files')),
  stream: (src) => (src.stream ? (providers()[src.stream] ?? null) : null),
};

/** Places to go: a gene's locus, an exon at letter zoom, a whole chromosome (its zoom levels and the gates' hints). */
const PLACES: [string, string, number, number][] = [
  ['TP53', 'chr17', 7_661_779, 7_687_538],
  ['TP53 exon 5', 'chr17', 7_675_990, 7_676_290],
  ['BRCA1', 'chr17', 43_044_295, 43_125_483],
  ['chr17 whole', 'chr17', 0, 83_257_441],
  ['CFTR', 'chr7', 117_480_025, 117_668_665],
  ['HBB', 'chr11', 5_225_464, 5_229_395],
];

function Panel({ chr, onChromosome, shimmerFps }: { chr: string; onChromosome(p: number): void; shimmerFps: number }) {
  const { context: ctx } = useSeqViewer();
  useLiveVersion(ctx.streams.stats);
  const st = ctx.streams.stats.current;
  useEffect(() => ctx.config.set(SeqViewerConfig.Stream.ShimmerFps, shimmerFps), [ctx, shimmerFps]);
  // On this chromosome: zoomed there; on another: its document (one chromosome a document, the navigator spanning it).
  const go = (k: number) => {
    const [, name, a, e] = PLACES[k];
    if (name === chr) ctx.viewports.get(JBROWSE_SECTION)?.zoomTo(0, a, e);
    else onChromosome(k);
  };
  const retry = () => void ctx.commands.run({ command: SeqCommands.Retry.id, params: {} }, 'api');
  return (
    <div className={f.stack}>
      <div className={f.row}>
        {PLACES.map(([label], k) => (
          <TextAction key={label} onClick={() => go(k)}>
            {label}
          </TextAction>
        ))}
        <TextAction onClick={retry}>Retry failed regions</TextAction>
      </div>
      <div className={f.code}>
        {st.requests} regions asked, {st.inFlight} in flight, {st.cancelled} cancelled, {st.errors} failed;{' '}
        {(st.bytes / 1024 / 1024).toFixed(1)} MiB held
      </div>
    </div>
  );
}

function Card({ shimmerFps }: { shimmerFps: number }) {
  const [place, setPlace] = useState(0);
  const [, chr, a, e] = PLACES[place];
  const d = useMemo(() => jbrowseDoc(chr, a, e), [chr, a, e]);
  const options = useMemo(
    () => ({
      loaders,
      maxHeight: 520,
      navigator: 'ideogram',
      registry: IDEOGRAM_REGISTRY,
      parts: ['TrackActions'] as const,
    }),
    [],
  );
  return (
    <ViewerCard
      key={chr}
      doc={d}
      options={options}
      label={`hg38 ${chr}: the reference, GENCODE, ClinVar and GC content`}
    >
      <Panel chr={chr} onChromosome={setPlace} shimmerFps={shimmerFps} />
    </ViewerCard>
  );
}

export const Live: Story = {
  name: 'hg38, live',
  args: { shimmerFps: 30 },
  argTypes: {
    shimmerFps: {
      control: { type: 'range', min: 0, max: 60, step: 5 },
      description: "`Stream.ShimmerFps`: the skeletons' shimmer, frames a second (0: a still band)",
    },
  },
  render: (args) => (
    <PrismPage
      title="The JBrowse 2 demo: hg38, read live from the public files"
      caption={
        <>
          Four indexed files read by byte range as you look, through `@molstar/seq-core/stream/gmod`: UCSC's hg38 2bit
          (the reference letters, asked only where a pixel holds 16 positions or fewer), GENCODE v50's BigBed
          (transcripts folded into a lane a gene while their lanes would not fit, exons joined, labeled by gene and
          colored by class; the track's … menu has Expand rows, a lane each), ClinVar's tabix VCF (variants colored and
          shaped by type; it names chromosomes `1`, `MT`, mapped to `chr1`, `chrM`) and UCSC's GC BigWig (its zoom
          levels drawn coarse under fine). Needs the network; nothing is cached but what the view holds.
        </>
      }
    >
      <div className={f.kicker}>
        It opens on TP53 (chr17), the navigator spanning the chromosome as its ideogram (hg38's cytobands, the window
        over them). Try its exon 5 at letter zoom; BRCA1 for ClinVar's dense regions (estimated counts first, from the
        index, "zoom in to see variants, or click to load them": a click loads what is in view; each variant once a
        pixel holds 16 positions or fewer); the whole chromosome, where GC and genes come from the files' zoom levels,
        ClinVar from its index's estimate, and the reference says to zoom in: nothing of the sequence and no item is
        fetched there. Hover a transcript or a variant for its fields; click one to select its range. Pending regions
        shimmer at the knob's rate; a failed region is a red band: click it to retry.
      </div>
      <Card shimmerFps={args.shimmerFps} />
    </PrismPage>
  ),
};
