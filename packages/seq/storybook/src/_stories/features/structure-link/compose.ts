/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { ChainTableFile as P53ChainTableFile, StructureAddresses } from '@molstar/seqx-examples/p53/structure';
import { ComplexAddress } from '@molstar/seqx-examples/complex/units';
import { readChainTable } from '@molstar/seqx-examples/p53/structure-record';
import { type LinkStructureOptions, linkStructure, type StructureBinding } from '@molstar/seqx-viewer/link-structure';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { ChainTable } from '../../support/fake-structure/chain-table.js';
import { FakeStructureViewer, type Selectors } from '../../support/fake-structure/fake-viewer.js';
import { colorByLayer, type HostAdapter, hostAdapter, viewerFacts } from '../../support/host-adapter.js';

// The structure link's stories, DOM-free (the stories and their headless twins call these): a fake structure viewer on
// a recorded chain table, its host's adapter over the viewer's complex, and `linkStructure` between them, with a log of
// the messages that crossed (what the story prints and a twin reads).

/** The datasets the fake has a chain table of: the entry's file (relative to the examples) and the complex's address. */
export const Datasets = {
  p53: { tableFile: P53ChainTableFile, complex: StructureAddresses.complex as string },
  complex: { tableFile: 'fixtures/complex/1FIN-chains.json', complex: ComplexAddress as string },
} as const;
export type DatasetName = keyof typeof Datasets;

/** The chain table of a dataset from its file's text (the story fetches it, a twin reads it from disk). */
export const chainTableOf = (text: string): ChainTable => readChainTable(text);

/** Selectors as a line: `A 96–105`, `E` (a ligand's chip), joined. */
export function describeSelectors(selectors: Selectors | null): string {
  if (!selectors?.length) return 'nothing';
  return selectors
    .map((s) => {
      const id = `${s.label_asym_id ?? '*'}${s.instance_id ? ` ${s.instance_id}` : ''}`;
      if (s.label_seq_id !== undefined) return `${id} ${s.label_seq_id}`;
      if (s.beg_label_seq_id !== undefined) return `${id} ${s.beg_label_seq_id}–${s.end_label_seq_id}`;
      return id;
    })
    .join(', ');
}

/** What Seq* lit in the fake, as residue numbers by chain (`A 161, B 161`): the readout where a one-cell shift is read. */
export function litLine(fake: FakeStructureViewer): string {
  const spans = fake.highlight.value;
  if (!spans.length) return 'Lit in the fake: nothing';
  const cells = spans.map((s) => {
    const row = fake.rows[s.row];
    const first = row.residues[s.start]?.auth ?? '?';
    return s.end - s.start > 1
      ? `${row.asym} ${first}–${row.residues[s.end - 1]?.auth ?? '?'}`
      : `${row.asym} ${first}`;
  });
  return `Lit in the fake: ${cells.join(', ')}`;
}

export interface FakeLinkOptions {
  /** Milliseconds the fake delays its messages (0). */
  readonly latency?: number;
  /** The link's options (its complexes, directions, limits). */
  readonly link?: LinkStructureOptions;
}

export interface FakeLink {
  readonly fake: FakeStructureViewer;
  readonly adapter: HostAdapter;
  /** The last messages in either direction, newest last: `→` what Seq* told the fake, `←` what the fake told. */
  readonly log: State<readonly string[]>;
  /** Adds a line of the host's own (the color command's). */
  say(line: string): void;
  /** Ends the link and the fake. */
  dispose(): void;
}

const LogLines = 10;

/**
 * The fake structure viewer of a chain table linked to a settled viewer's complex. The adapter is made from the
 * complex's facts (an address that is no complex throws), the side logs both ways, the link clears with the fake.
 */
export function linkFake(
  viewer: SeqViewer,
  table: ChainTable,
  complex: string,
  options: FakeLinkOptions = {},
): FakeLink {
  const facts = viewerFacts(viewer, complex);
  if (!facts) throw new Error(`${complex} is no complex of the view`);
  const fake = new FakeStructureViewer(table, { latency: options.latency ?? 0 });
  const adapter = hostAdapter(fake, facts);
  const log = new State<readonly string[]>([]);
  const push = (line: string) => log.set([...log.value.slice(1 - LogLines), line]);
  const { side } = adapter;
  const logged: typeof side = {
    ...side,
    highlight: (t) => {
      push(`→ highlight ${describeSelectors(t && adapter.selectorsOfTarget(t))}`);
      side.highlight(t);
    },
    select: (t) => {
      push(`→ select ${describeSelectors(t && adapter.selectorsOfTarget(t))}`);
      side.select(t);
    },
    focus: (t, o) => {
      push(`→ focus ${describeSelectors(adapter.selectorsOfTarget(t))} over ${o.durationMs} ms`);
      side.focus(t, o);
    },
  };
  const told = (['hover', 'select', 'focus'] as const).map((kind) =>
    fake[kind === 'hover' ? 'onHover' : kind === 'select' ? 'onSelect' : 'onFocus']((s) =>
      push(`← ${kind} ${describeSelectors(s)}`),
    ),
  );
  const unlink = linkStructure(viewer, logged, options.link);
  return {
    fake,
    adapter,
    log,
    say: push,
    dispose() {
      unlink();
      for (const off of told) off();
      fake.dispose();
    },
  };
}

/**
 * What the host's color command runs (`colorByBehavior`): the viewer's track colored onto the current link's fake, the
 * message logged as a line of the link's own (`→ color …`).
 */
export function colorRunner(
  viewer: SeqViewer,
  current: () => FakeLink | undefined,
  scheme: () => 'light' | 'dark',
  bindings: () => readonly StructureBinding[] | undefined = () => undefined,
): (layer: string) => string {
  return (layer) => {
    const l = current();
    if (!l) return 'The structure is not linked yet';
    const message = colorByLayer(viewer, l.fake, l.adapter, layer, scheme(), bindings());
    l.say(`→ color ${layer}: ${message}`);
    return message;
  };
}
