/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { largeAlignment, longSequence } from "./support/seqstar-docs-views.js";
import { GenomeView } from "./support/seqstar-genome-view.js";
import { PrismPage } from "./support/seqstar-prism.js";
import { ViewerCard } from "./support/seqstar-viewer.js";
import { type Args, viewerControls } from "./support/seqstar-viewer-args.js";

/**
 * Large documents: alignments and sequences too big to read at once, drawn from reductions, and a genome of 3,000,000
 * nucleotides.
 */
const meta = {
  title: "Seq*/Performance/Large documents",
  ...viewerControls,
  // Said here too: Storybook reads a story file's layout from its own meta, not through a spread.
  parameters: { ...viewerControls.parameters, layout: "fullscreen" },
} satisfies Meta<Args>;
export default meta;
type Story = StoryObj<Args>;

const LARGE = new Map<string, SeqViewSpecDoc>();
const large = (rows: number, columns: number) => {
  const key = `${rows}x${columns}`;
  const d = LARGE.get(key) ?? largeAlignment(rows, columns);
  LARGE.set(key, d);
  return d;
};
function Large({ rows, columns, a }: { rows: number; columns: number; a: Args }) {
  return (
    <PrismPage
      title={
        <>
          A synthetic alignment of {rows.toLocaleString("en-US")} rows × {columns.toLocaleString("en-US")} columns
        </>
      }
      caption={
        <>
          Synthetic, by value. Zoom out (Ctrl-wheel) and Alt-wheel rows: below a device pixel the rows draw from pyramid
          reductions (each sample the majority class, opaque by its share), computed in place by the tile manager,
          coarse under fine.
        </>
      }
    >
      <ViewerCard
        doc={large(rows, columns)}
        label="Alignment"
        options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }}
      />
    </PrismPage>
  );
}
/** 10,000 rows × 300 columns in memory: reductions computed by the tile manager when zoomed out. */
export const LargeAlignment: Story = {
  name: "Large alignment (10,000 × 300)",
  render: (a) => <Large rows={10_000} columns={300} a={a} />,
};
/** 1,000 rows × 2,000 columns. */
export const WideAlignment: Story = {
  name: "Large alignment (1,000 × 2,000)",
  render: (a) => <Large rows={1_000} columns={2_000} a={a} />,
};

let LONG: SeqViewSpecDoc | null = null;
export const LongSequence: Story = {
  name: "35,000 residues",
  render: (a) => {
    LONG ??= longSequence(35_000);
    return (
      <PrismPage
        title="A titin-sized sequence"
        caption="Synthetic, by value: 35,000 residues and 70 domains on one axis."
      >
        <ViewerCard doc={LONG} label="Sequence and domains" options={{ navigator: a.navigator, ui: a.ui }} />
      </PrismPage>
    );
  },
};

/**
 * A genome at scale: 3,000,000 nucleotides in 24 chromosome blocks, genes as
 * transcripts, proteins over their codons, GC content, variants and repeats; synthetic, seeded, made in the page.
 */
export const GenomeScale: Story = {
  name: "Genome-scale, 3M nucleotides",
  render: (a) => <GenomeView navigator={a.navigator} ui={a.ui} />,
};
