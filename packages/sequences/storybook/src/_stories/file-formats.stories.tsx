/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryTable, MolstarSvsEntry, readStructure, type StructureFormat } from '@molstar/seq-molstar-poc';
import { useLiveVersion, useSeqState, useSeqViewer } from '@molstar/seq-viewer/react';
import type { Meta, StoryObj } from '@storybook/react-vite';
import { type ReactNode, Suspense, use, useEffect, useMemo, useState } from 'react';
import kinaseSvsx from '../../../core/_test/format/fixtures/kinase.svsx?url';
import featuresJson from '../../../molstar-poc/_test/fixtures/p53/P04637.json?raw';
import variantsJson from '../../../molstar-poc/_test/fixtures/p53/variation-P04637.json?raw';
import ocjUrl from '../../../molstar-poc/_test/fixtures/structures/2ocj.bcif?url';
import afUrl from '../../../molstar-poc/_test/fixtures/structures/AF-P04637-F1-model_v6.bcif?url';
import ins1 from '../../../molstar-poc/_test/fixtures/structures/ins1.cif?raw';
import {
  AF_P53_URI,
  ALIGNMENT_FORMATS,
  type AlignmentFormat,
  alignmentCase,
  archiveCase,
  DOCUMENT_FORMATS,
  type DocumentFormat,
  documentCase,
  type Fetched,
  type FormatCase,
  GENOME_HOST_CODE,
  type Given,
  genomeCase,
  liveLoaders,
  nodesText,
  P53_URIS,
  PDBE_AFDB_SHOWN,
  pdbeAfdbDoc,
  pdbeAfdbFiles,
  recordedLoaders,
  STRUCTURE_FORMATS,
  type StructureForm,
  sequencesCase,
  structuresCase,
  TABLE_FORMATS,
  type TableFormat,
  TILES_URI,
  tablesCase,
  tilesCase,
  UNIPROT_SHOWN,
  uniprotDoc,
  uniprotFiles,
} from './support/seqstar-file-formats.js';
import f from './support/seqstar-format.module.css';
import { IssueList } from './support/seqstar-issues.js';
import { PrismPage, prism } from './support/seqstar-prism.js';
import { ViewerCard, type ViewerOptions } from './support/seqstar-viewer.js';
import { Choices, ScrollArea } from './ui/index.js';

/**
 * Loading data, one story a family of formats (seqspec:sources#formats-and-services): the nodes that load it as the
 * document holds them, the data as given (a file's start, or the requests its loaders answered), the viewer drawing
 * it, and what loading reported. Offline by default, on the committed fixtures or files written here from them; the
 * service stories have a live mode that fetches the public APIs. Documents and loaders: seqstar-file-formats.ts.
 */
const meta = { title: 'Seq*/File formats', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj<typeof meta>;

/** The Mol* adapter's format half, for the stories that read structures (its formats and services, real parsers). */
const MOLSTAR = [{ svs: MolstarSvsEntry }];

/** A code or data block with its title. */
function Code({ label, text }: { label: ReactNode; text: string }) {
  return (
    <div className={f.stack}>
      <h3 className={f.subheading}>{label}</h3>
      <ScrollArea className={f.box} orientation="both">
        <pre className={f.pre}>{text}</pre>
      </ScrollArea>
    </div>
  );
}

/** What the viewer reported while loading (notes included), under it. */
function Issues() {
  const { context } = useSeqViewer();
  const issues = useSeqState(context.issues.current);
  return (
    <div className={f.stack}>
      <h3 className={f.subheading}>What loading reported</h3>
      <IssueList issues={issues} />
    </div>
  );
}

/** The requests the loaders answered, read again as the graph settles. */
function Requests({ log }: { log: readonly Fetched[] }) {
  const { context } = useSeqViewer();
  useSeqState(context.graph.assembled);
  useSeqState(context.issues.current);
  return (
    <div className={f.stack}>
      <h3 className={f.subheading}>The requests, as answered</h3>
      {log.length === 0 ? <div className={f.none}>none yet</div> : null}
      {log.map((x) => (
        <Code key={x.uri} label={`${x.uri}, ${x.bytes.toLocaleString('en-US')} bytes`} text={x.head} />
      ))}
    </div>
  );
}

/** The document's nodes that load the data, and the data as given. */
function Sources({ c, given = c.given, children }: { c: FormatCase; given?: readonly Given[]; children?: ReactNode }) {
  const nodes = useMemo(() => nodesText(c.doc, c.shown), [c]);
  return (
    <div className={prism.pair}>
      <section className={f.card}>
        <h2 className={f.heading}>In the document</h2>
        <Code label="The nodes that load it (long lists abridged)" text={nodes} />
      </section>
      <section className={f.card}>
        <h2 className={f.heading}>The data as given</h2>
        {given.map((g) => (
          <Code key={g.label} label={g.label} text={g.text} />
        ))}
        {children}
      </section>
    </div>
  );
}

/** A row of choices in a card of its own, above the viewer. */
function Pick<T extends string>({
  label,
  choices,
  active,
  onPick,
}: {
  label: string;
  choices: readonly { id: T; label: string }[];
  active: T;
  onPick(id: T): void;
}) {
  return (
    <section className={f.card}>
      <div className={prism.controls}>
        <div className={prism.wide}>
          <Choices label={label} choices={choices} active={active} onPick={(id) => onPick(id as T)} />
        </div>
      </div>
    </section>
  );
}

/** The case's viewer with what loading reported, and `children` (the requests) under it. */
function Viewer({
  c,
  label,
  options,
  children,
}: {
  c: FormatCase;
  label: ReactNode;
  options?: ViewerOptions;
  children?: ReactNode;
}) {
  const o = useMemo(() => ({ loaders: c.loaders, maxHeight: 420, ...options }), [c, options]);
  return (
    <ViewerCard doc={c.doc} label={label} options={o}>
      <div className={prism.actions}>
        <Issues />
        {children}
      </div>
    </ViewerCard>
  );
}

// ---- Sequences ----

function SequencesPage() {
  const c = useMemo(sequencesCase, []);
  return (
    <PrismPage
      title="Sequences: FASTA and residues"
      caption={
        <>
          Three ways a sequence reads: from a FASTA file by URI, one record each by its id (<code>record</code>); from a
          FASTA in place (a source's <code>value</code>: the whole file, its header kept); and the residues themselves
          (a sequence's <code>value</code>). The file's last record has a header and no sequence: it is read, with a
          note (<code>svs.source.empty_record</code>). Barnase and barstar from UniProt's P00648 and P11540 (CC BY 4.0),
          the construct and the design are the format's test fixtures. Recorded, no network.
        </>
      }
    >
      <Viewer c={c} label="Four sequences, four blocks: the alphabet detected (the construct is DNA)" />
      <Sources c={c} />
    </PrismPage>
  );
}
export const Sequences: Story = { name: 'Sequences', render: () => <SequencesPage /> };

// ---- Alignments and trees ----

const ALIGNMENT_LABELS: Record<AlignmentFormat, string> = {
  fasta: 'Aligned FASTA',
  a3m: 'A3M',
  stockholm: 'Stockholm',
  clustal: 'Clustal',
};

function AlignmentsPage() {
  const [format, setFormat] = useState<AlignmentFormat>('stockholm');
  const c = useMemo(() => alignmentCase(format), [format]);
  return (
    <PrismPage
      title="Alignments and trees: aligned FASTA, A3M, Stockholm, Clustal, Newick"
      caption={
        <>
          Pfam's PF00042 seed (its first eight globins, CC0) in four alignment formats, and a Newick tree over them
          (written for the story), matched to the rows by name: a leaf <code>HBA_HUMAN</code> is the row{' '}
          <code>HBA_HUMAN/27-137</code>. Aligned FASTA reads every column; A3M reads lower case as inserts, outside the
          columns (the first row's gaps made into inserts here; noted as <code>svs.alignment.inserts</code>); Stockholm
          keeps its tags, and its <code>#=GS … AC</code> lines are read as a rows table (<code>tag: "GS:AC"</code>, a
          data column). The Stockholm file was cut before its closing <code>{'//'}</code>: it reads, with a note (
          <code>svs.source.unterminated</code>). Clustal (<code>.aln</code>, CLUSTAL W's layout) joins each row's blocks
          in order; its residue counts and conservation lines are not read.
        </>
      }
    >
      <Pick<AlignmentFormat>
        label="The alignment's format"
        choices={ALIGNMENT_FORMATS.map((id) => ({ id, label: ALIGNMENT_LABELS[id] }))}
        active={format}
        onPick={setFormat}
      />
      <Viewer key={format} c={c} label={`PF00042, eight rows, from ${ALIGNMENT_LABELS[format]}`} />
      <Sources c={c} />
    </PrismPage>
  );
}
export const Alignments: Story = { name: 'Alignments and trees', render: () => <AlignmentsPage /> };

// ---- Tables ----

const TABLE_LABELS: Record<TableFormat, string> = {
  tsv: 'TSV',
  csv: 'CSV',
  json_rows: 'JSON, rows',
  json_columns: 'JSON, columns',
};

function TablesPage() {
  const [format, setFormat] = useState<TableFormat>('tsv');
  const c = useMemo(() => tablesCase(format), [format]);
  return (
    <PrismPage
      title="Tables: TSV, CSV and JSON"
      caption={
        <>
          Two annotations read from table files, in four forms: TSV and CSV (cells are text, typed by their column's
          role), and JSON as an array of row objects or an object of columns (both read the same). The files name their
          own columns; <code>fields</code> maps them to the format's (<code>first</code> → <code>beg</code>,{' '}
          <code>kd</code> → <code>value</code>). Positions are Kabat labels: the hydropathy file gives each residue's
          number and insertion code apart, and <code>ins_code</code> joins them (52 and A name 52A), so the inserted
          residues 52A, 82A–C and 100A–C get their own values. Trastuzumab's VH from the format's fixtures; the files
          are written for the story. Recorded, no network.
        </>
      }
    >
      <Pick<TableFormat>
        label="The tables' format"
        choices={TABLE_FORMATS.map((id) => ({ id, label: TABLE_LABELS[id] }))}
        active={format}
        onPick={setFormat}
      />
      <Viewer key={format} c={c} label={`Trastuzumab VH in Kabat numbering, its tables from ${TABLE_LABELS[format]}`} />
      <Sources c={c} />
    </PrismPage>
  );
}
export const Tables: Story = { name: 'Tables', render: () => <TablesPage /> };

// ---- Structure files ----

/** A CIF category's first rows as a text table, from a structure file read with the adapter's parser. */
function useCategory(bytes: Uint8Array | null, format: StructureFormat, name: string, rows = 10): string {
  const [text, setText] = useState('reading…');
  useEffect(() => {
    if (!bytes) return;
    let live = true;
    void readStructure(bytes, format).then((src) => {
      const t = src.frame ? categoryTable(src.frame, name) : null;
      if (!live) return;
      if (!t) return setText(`no category ${name}`);
      const names = Object.keys(t.columns);
      const cells = [
        names,
        ...Array.from({ length: Math.min(rows, t.length) }, (_, i) => names.map((k) => t.columns[k][i] ?? '?')),
      ];
      const w = names.map((_, k) => Math.max(...cells.map((r) => r[k].length)));
      const lines = cells.map((r) =>
        r
          .map((x, k) => x.padEnd(w[k]))
          .join('  ')
          .trimEnd(),
      );
      setText(`${lines.join('\n')}\n… (${t.length.toLocaleString('en-US')} rows)`);
    });
    return () => {
      live = false;
    };
  }, [bytes, format, name, rows]);
  return text;
}
/** A fixture's bytes by its bundled URL, once. */
function useBytes(url: string): Uint8Array | null {
  const [bytes, setBytes] = useState<Uint8Array | null>(null);
  useEffect(() => {
    let live = true;
    void fetch(url)
      .then((r) => r.arrayBuffer())
      .then((b) => live && setBytes(new Uint8Array(b)));
    return () => {
      live = false;
    };
  }, [url]);
  return bytes;
}

const STRUCTURE_LABELS: Record<StructureForm, string> = {
  bcif: 'BinaryCIF: an AlphaFold DB model',
  mmcif: 'mmCIF: insertion codes',
};
const STRUCTURE_OPTIONS: ViewerOptions = { registry: MOLSTAR, fixtureStructures: false };

function StructureCase({ form, af }: { form: StructureForm; af: Uint8Array }) {
  const c = useMemo(() => structuresCase(form, { af, ins1 }), [form, af]);
  const category = useCategory(form === 'bcif' ? af : null, 'bcif', 'ma_qa_metric_local');
  const given: readonly Given[] =
    form === 'bcif'
      ? [{ label: `${AF_P53_URI} (recorded): its ma_qa_metric_local, as Mol*'s parser reads it`, text: category }]
      : c.given;
  return (
    <>
      <Viewer
        key={form}
        c={c}
        label={form === 'bcif' ? 'AF-P04637-F1, chain A' : 'INS1, chain H: 52, 52A, 52B in auth numbering'}
        options={STRUCTURE_OPTIONS}
      />
      <Sources c={c} given={given} />
    </>
  );
}

function StructuresPage() {
  const [form, setForm] = useState<StructureForm>('bcif');
  const af = useBytes(afUrl);
  return (
    <PrismPage
      title="Structure files: mmCIF and BinaryCIF"
      caption={
        <>
          A structure file read by URI through the Mol* adapter (<code>@molstar/seq-molstar-poc</code>, Mol*'s own
          parsers, headless): the <code>molstar:structure</code> complex provider makes its chains units, and a CIF
          category is read as a table (<code>from_source</code>, <code>category_name</code>), its columns mapped by{' '}
          <code>fields</code>. AlphaFold DB's model of P04637 (BinaryCIF, CC BY 4.0) gives its pLDDT from{' '}
          <code>ma_qa_metric_local</code> by <code>label_seq_id</code>; a small hand-written mmCIF gives B-factors from{' '}
          <code>atom_site</code> by author numbers, <code>auth_seq_id</code> with <code>pdbx_PDB_ins_code</code> as{' '}
          <code>ins_code</code>, so 52A and 52B are their own residues. PDB, SDF and MOL2 read the same way, without
          categories. Recorded, no network.
        </>
      }
    >
      <Pick<StructureForm>
        label="The file"
        choices={STRUCTURE_FORMATS.map((id) => ({ id, label: STRUCTURE_LABELS[id] }))}
        active={form}
        onPick={setForm}
      />
      {af ? <StructureCase form={form} af={af} /> : <div className={f.status}>Reading the fixture…</div>}
    </PrismPage>
  );
}
export const Structures: Story = { name: 'Structure files', render: () => <StructuresPage /> };

// ---- Services ----

type Mode = 'recorded' | 'live';
const MODES: readonly { id: Mode; label: string }[] = [
  { id: 'recorded', label: 'Recorded' },
  { id: 'live', label: 'Live (calls the public API)' },
];

function UniprotPage() {
  const [mode, setMode] = useState<Mode>('recorded');
  const doc = useMemo(uniprotDoc, []);
  const c = useMemo<FormatCase>(
    () => ({
      doc,
      loaders: mode === 'live' ? liveLoaders() : recordedLoaders(uniprotFiles(featuresJson, variantsJson)),
      shown: UNIPROT_SHOWN,
      given: [],
    }),
    [doc, mode],
  );
  return (
    <PrismPage
      title="UniProt: sequence, features, variants"
      caption={
        <>
          The kernel's one service, <code>uniprot</code>: a record id and a resource resolve to a request,{' '}
          <code>sequence</code> to UniProt's FASTA (the sequence gets the numbering <code>uniprot:P04637</code>),{' '}
          <code>features</code> to the entry's JSON and <code>variants</code> to the Proteins API's variation, both read
          as ranges. Recorded: the responses committed for P04637 (UniProtKB, CC BY 4.0; trimmed to what the readers
          take). Live calls the public APIs from your browser (rest.uniprot.org, www.ebi.ac.uk), the whole entry and its
          thousands of variants; those starting past the sequence's end (a stop lost, a frameshift's new tail) are left
          out and those running past it clipped, one note (<code>svs.source.skipped</code>) counting them.
        </>
      }
    >
      <Pick<Mode> label="Source" choices={MODES} active={mode} onPick={setMode} />
      <Viewer key={mode} c={c} label={`P04637 (p53), ${mode}`}>
        <Requests log={c.loaders.fetched} />
      </Viewer>
      <Sources
        c={c}
        given={[{ label: "The service's requests", text: Object.values(P53_URIS).slice(0, 3).join('\n') }]}
      />
    </PrismPage>
  );
}
export const UniProt: Story = { name: 'UniProt', render: () => <UniprotPage /> };

const PDBE_OPTIONS: ViewerOptions = { registry: MOLSTAR, fixtureStructures: false };

function PdbeAfdbCase({ mode, pdbe, afdb }: { mode: Mode; pdbe: Uint8Array; afdb: Uint8Array }) {
  const doc = useMemo(pdbeAfdbDoc, []);
  const c = useMemo<FormatCase>(
    () => ({
      doc,
      loaders: mode === 'live' ? liveLoaders() : recordedLoaders(pdbeAfdbFiles({ pdbe, afdb })),
      shown: PDBE_AFDB_SHOWN,
      given: [],
    }),
    [doc, mode, pdbe, afdb],
  );
  const sifts = useCategory(pdbe, 'bcif', 'pdbx_sifts_unp_segments');
  return (
    <>
      <Viewer key={mode} c={c} label={`AF-P04637-F1 with 2OCJ's chains, ${mode}`} options={PDBE_OPTIONS}>
        <Requests log={c.loaders.fetched} />
      </Viewer>
      <Sources
        c={c}
        given={[
          { label: "The services' requests", text: Object.values(P53_URIS).slice(3).join('\n') },
          {
            label: '2ocj.bcif (recorded): its pdbx_sifts_unp_segments, the SIFTS the numbering comes from',
            text: sifts,
          },
        ]}
      />
    </>
  );
}

function PdbeAfdbPage() {
  const [mode, setMode] = useState<Mode>('recorded');
  const pdbe = useBytes(ocjUrl);
  const afdb = useBytes(afUrl);
  return (
    <PrismPage
      title="PDBe and AlphaFold DB: structures, SIFTS, pLDDT and PAE"
      caption={
        <>
          The Mol* adapter's services: <code>pdbe</code> (<code>structure</code>, PDBe's BinaryCIF with the updated
          mmCIF's SIFTS; <code>sifts</code>, the updated mmCIF) and <code>afdb</code> (<code>structure</code>, the
          model's BinaryCIF; <code>pae</code>, its predicted aligned error). 2OCJ's chains carry the numbering{' '}
          <code>uniprot:P04637</code> from SIFTS and map onto P04637 by it; AlphaFold DB's model maps onto P04637 one to
          one, its pLDDT from its B-factors (<code>molstar:plddt</code>), its PAE a dot plot (<code>molstar:pae</code>,
          in the closed group), 2OCJ's coverage placed on it through P04637. Recorded: PDBe's and AlphaFold DB's files
          committed (wwPDB CC0, AlphaFold DB CC BY 4.0) and a <em>synthetic</em> PAE made in the page (AlphaFold DB's is
          about 400 KB). Live calls the public files from your browser (www.ebi.ac.uk, alphafold.ebi.ac.uk,
          rest.uniprot.org).
        </>
      }
    >
      <Pick<Mode> label="Source" choices={MODES} active={mode} onPick={setMode} />
      {pdbe && afdb ? (
        <PdbeAfdbCase mode={mode} pdbe={pdbe} afdb={afdb} />
      ) : (
        <div className={f.status}>Reading the fixtures…</div>
      )}
    </PrismPage>
  );
}
export const PdbeAfdb: Story = { name: 'PDBe and AlphaFold DB', render: () => <PdbeAfdbPage /> };

// ---- Streamed genome files ----

/** Each file's reads so far: how many, their bytes, the last few ranges. */
function Reads({ reads }: { reads: Readonly<Record<string, [number, number][]>> }) {
  const { context } = useSeqViewer();
  useLiveVersion(context.streams.stats);
  const st = context.streams.stats.current;
  const lines = Object.entries(reads).map(([file, list]) => {
    const bytes = list.reduce((n, [, l]) => n + l, 0);
    const last = list.slice(-3).map(([o, l]) => `${o}–${o + l - 1}`);
    return `${file.padEnd(7)} ${String(list.length).padStart(3)} reads, ${bytes.toLocaleString('en-US').padStart(7)} bytes; last: ${last.join(', ') || 'none'}`;
  });
  return (
    <Code
      label={`The byte ranges read: ${st.requests} regions asked of the providers, ${st.inFlight} in flight, ${st.errors} failed`}
      text={lines.join('\n')}
    />
  );
}

function GenomePage() {
  const c = useMemo(genomeCase, []);
  return (
    <PrismPage
      title="Streamed genome files: 2bit, BigWig, BigBed, tabix VCF"
      caption={
        <>
          Indexed genome files read by region through the host's stream providers (
          <code>@molstar/seq-core/stream/gmod</code>), never whole: the sources name a stream and its{' '}
          <code>data_kind</code>, the host maps each to a reader over byte ranges. A synthetic genome of two sequences (
          <code>chrS1</code>, <code>chrS2</code>), its files written in memory by the format's test helper (
          <code>demo-files.ts</code>): a 2bit (an N run, soft-masked repeats), a BigWig of GC percent, a BigBed of BED12
          genes (exons joined by id) and a bgzipped VCF with its tabix index. Pan and zoom to see the reads grow; zoom
          in for letters. hg38's real files, live, are the{' '}
          <a href="?path=/story/seq-cases-genome-browser--live">genome browser demo</a>.
        </>
      }
    >
      <Viewer c={c} label="chrS1 of the synthetic genome, by region">
        <Reads reads={c.reads} />
      </Viewer>
      <Sources c={c}>
        <Code label="The host: a provider per file (as in this story)" text={GENOME_HOST_CODE} />
      </Sources>
    </PrismPage>
  );
}
export const GenomeFiles: Story = { name: 'Streamed genome files', render: () => <GenomePage /> };

// ---- The tiles server ----

/** The tile server's requests so far: how many of each resource, and the first paths. */
function Calls({ calls }: { calls: readonly string[] }) {
  const { context } = useSeqViewer();
  useLiveVersion(context.streams.stats);
  const kinds = new Map<string, number>();
  for (const p of calls) kinds.set(p.split('/')[0], (kinds.get(p.split('/')[0]) ?? 0) + 1);
  const summary = [...kinds].map(([k, n]) => `${k} ${n}`).join(', ');
  return (
    <Code
      label={`${calls.length} requests to ${TILES_URI}${summary ? `: ${summary}` : ''}`}
      text={calls.slice(0, 24).join('\n') + (calls.length > 24 ? `\n… ${calls.length - 24} more` : '')}
    />
  );
}

function TilesPage() {
  const c = useMemo(tilesCase, []);
  return (
    <PrismPage
      title="The tiles protocol: seqstar-tiles/1"
      caption={
        <>
          A source with <code>protocol: "seqstar-tiles/1"</code> and a <code>uri</code> is streamed by the tiles client
          (<code>@molstar/seq-core/stream/tiles</code>): an alignment's 256 × 256 samples at a level of rows and columns
          (<code>tiles/ly/lx/ty/tx</code>, each sample the most frequent residue and its share), row labels (
          <code>labels/a/b</code>), whole rows (<code>rows/k</code>) and a column statistic (<code>stats/…</code>, here
          the conservation <code>column_stats</code> asks for). The document declares the alignment's shape; the
          manifest is built from it. Served by the in-process reference server (<code>tileServer</code>, 4,000
          procedural rows): the host gives the client its transport, HTTP GET under the URI by default.
        </>
      }
    >
      <Viewer c={c} label="4,000 rows × 360 columns, by tiles">
        <Calls calls={c.calls} />
      </Viewer>
      <Sources c={c} />
    </PrismPage>
  );
}
export const Tiles: Story = { name: 'The tiles server', render: () => <TilesPage /> };

// ---- Seq*'s own documents ----

const DOCUMENT_LABELS: Record<DocumentFormat, string> = {
  svsj: '.svsj (JSON)',
  svsm: '.svsm (MessagePack)',
  svsx: '.svsx (an archive)',
};

/** The fixture archive, read once. */
let archive: Promise<Awaited<ReturnType<typeof archiveCase>>> | null = null;
const readArchive = () =>
  (archive ??= fetch(kinaseSvsx)
    .then((r) => r.arrayBuffer())
    .then((b) => archiveCase(new Uint8Array(b))));

function ArchiveCase() {
  const { c, svsx } = use(readArchive());
  return (
    <>
      <Viewer c={c} label={`The document read from ${DOCUMENT_LABELS.svsx}, over its own files`}>
        <div className={f.kicker}>
          Its sources name the files by paths relative to the archive's base (its main file's directory:{' '}
          {svsx.base || 'the root'}); the viewer's loaders, <code>archiveLoaders</code>, carry that base, so the
          document validates and each path is answered from the archive.
        </div>
        <Requests log={c.loaders.fetched} />
      </Viewer>
      <Sources c={c} />
    </>
  );
}

function DocumentsPage() {
  const [format, setFormat] = useState<DocumentFormat>('svsj');
  const c = useMemo(() => (format === 'svsx' ? null : documentCase(format)), [format]);
  return (
    <PrismPage
      title="Seq*'s own documents: .svsj, .svsm, .svsx"
      caption={
        <>
          A SeqViewSpec document as JSON (<code>.svsj</code>) or MessagePack (<code>.svsm</code>, typed arrays as its
          typed-array extension), told apart by their first byte (<code>readDocument</code>), and the <code>.svsx</code>{' '}
          archive: a zip with the document as <code>index.svsm</code> (or <code>index.svsj</code>) and the files it
          reads, by relative paths. One document in all three: an ABL kinase excerpt from a FASTA and its sites from a
          TSV (the format's <code>kinase.svsx</code> fixture); the JSON and MessagePack ones name the files by their
          absolute URIs, answered from recorded text. No network.
        </>
      }
    >
      <Pick<DocumentFormat>
        label="The document's format"
        choices={[...DOCUMENT_FORMATS, 'svsx' as const].map((id) => ({ id, label: DOCUMENT_LABELS[id] }))}
        active={format}
        onPick={setFormat}
      />
      {c ? (
        <>
          <Viewer key={format} c={c} label={`The document read from ${DOCUMENT_LABELS[format]}`} />
          <Sources c={c} />
        </>
      ) : (
        <Suspense fallback={<div className={f.status}>Opening the archive…</div>}>
          <ArchiveCase />
        </Suspense>
      )}
    </PrismPage>
  );
}
export const Documents: Story = { name: 'Seq* documents', render: () => <DocumentsPage /> };
