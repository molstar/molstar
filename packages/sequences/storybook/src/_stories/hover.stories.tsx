/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import type { Meta, StoryObj } from "@storybook/react-vite";
import type { ReactNode } from "react";
import { splitUnitsDoc } from "../../../viewer/_test/split-units-fixtures.ts";
import { logoDoc, tracksDoc } from "./support/seqstar-docs-views";
import f from "./support/seqstar-format.module.css";
import { PrismPage } from "./support/seqstar-prism";
import { ViewerCard } from "./support/seqstar-viewer";

/**
 * Hover emphasis: the element under the pointer, and only it, outlined inside its own pixels in an
 * ink derived from its color, darker on a light skin and lighter on a dark one. Ranges (a transcript's every exon),
 * points, stacked variants, secondary structure, values and a matrix's cell, in the toolbar's theme (a logo gets none: a frame
 * over its stretched letters read as wrong, and the hover band marks the column); the
 * keyboard cursor's item too (click the viewer, then the arrow keys); and a section with `hover_emphasis: false`.
 */
const meta = { title: "Seq*/Interaction/Hover", parameters: { layout: "fullscreen" } } satisfies Meta;
export default meta;
type Story = StoryObj;

/** One viewer, in the toolbar's theme (light, dark, blue). */
function One({ doc, title }: { doc: SeqViewSpecDoc; title: string }) {
  return <ViewerCard doc={doc} label={title} />;
}

function Frame({ title, kicker, children }: { title: string; kicker: ReactNode; children: ReactNode }) {
  return (
    <PrismPage title={title} caption={kicker}>
      <div className={f.stack}>{children}</div>
    </PrismPage>
  );
}

const TRACKS = tracksDoc("p53, every emphasized type");
const LOGO = logoDoc();
const OFF = tracksDoc("p53, hover emphasis off", { hover_emphasis: false });

export const Emphasis: Story = {
  name: "Hover emphasis",
  render: () => (
    <Frame
      title="Hover emphasis"
      kicker={
        <>
          Hover a range, a site, a variant, a strand, a score bar or a contact cell: the element under the pointer, and
          only it, is outlined inside its own pixels in an ink derived from its color (darker on light, lighter on dark;
          3:1 against the surface). A region in two parts (TAD) outlines both; a narrow item is filled in the ink; a
          glyph is ringed only where the ring fits its cell. The keyboard cursor's item too: click a viewer, then use
          the arrow keys. The fingerprint's stripes: <a href="?path=/story/seq-cases-viewer--case-6">case 6</a>.
        </>
      }
    >
      <One doc={TRACKS} title="Tracks" />
      <One doc={LOGO} title="A profile's cells (its logo: the band alone)" />
    </Frame>
  ),
};

export const Off: Story = {
  name: "Hover emphasis off",
  render: () => (
    <Frame
      title="Hover emphasis off"
      kicker={
        <>
          <span className={f.code}>
            "style": {"{"}"hover_emphasis": false{"}"}
          </span>{" "}
          on the section (or the render prop <span className={f.code}>hoverEmphasis</span>): the hover band alone, with
          no outline.
        </>
      }
    >
      <One doc={OFF} title="Off" />
    </Frame>
  ),
};

const SPLIT = splitUnitsDoc({ highlights: true });

export const RowMarks: Story = {
  name: "Per-row marks",
  render: () => (
    <Frame
      title="Per-row marks"
      kicker={
        <>
          A split track whose rows sit on units of their own: three copies of a chain (A, B, C), each observed in part,
          drawn on the reference's axis through a mapping, a row each. A mark on one copy lights its own row alone (and
          the reference's tracks): here a highlight on chain B's residues 4–7. A mark on the reference lights every row:
          the accent at 30–33. Hover a row: off its bars the position is that copy's residue, on a bar the bar's, and
          only that row lights; hover the reference's letters and every row lights. Rows split by any other field (one
          unit for all of them, or rows on the block's own sequence) light together, as before.
        </>
      }
    >
      <One doc={SPLIT} title="Chains on the reference" />
    </Frame>
  ),
};
