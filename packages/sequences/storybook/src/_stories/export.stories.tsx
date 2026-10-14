/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from "@storybook/react-vite";
import { Suspense, use, useMemo } from "react";
import { genomeDoc, synthGenome } from "../../../core/_test/format/genome-fixture.js";
import { compiledCase } from "./support/seqstar-cases.js";
import { ExportCard, ExportFrame } from "./support/seqstar-export.js";
import f from "./support/seqstar-format.module.css";

/**
 * Exporting a figure: the viewer's `view.exportSvg` action, which the UI's Export part runs, on
 * two live viewers: case 14 (a protein over its coding DNA) and the genome story (3,000,000 nucleotides). The SVG it makes
 * is drawn beside the canvas; the PNG is that SVG drawn by the page. File sizes, export times and downloads are shown.
 */
const meta = {
  title: "Seq*/Cases/Default UIs",
  parameters: { layout: "fullscreen", docs: { source: { type: "code" } } },
} satisfies Meta;
export default meta;
type Story = StoryObj<typeof meta>;

function Case14() {
  const c = use(compiledCase("case14"));
  return (
    <ExportCard
      doc={c.doc}
      label="Case 14: A protein over its coding DNA"
      note="A gene fragment of about 120 nucleotides, the protein translated from its codons, exons, a domain and a phosphosite."
    />
  );
}

/** The genome story's document, made once (synthetic, seeded; seconds to build). */
let GENOME: ReturnType<typeof genomeDoc> | null = null;
function Genome() {
  const doc = useMemo(() => {
    GENOME ??= genomeDoc(synthGenome());
    return GENOME;
  }, []);
  return (
    <ExportCard
      doc={doc}
      label="Genome-scale, 3M nucleotides"
      note="24 chromosomes, 3,000 genes, 300 proteins, GC content, 20,000 variants, 10,000 repeats: the figure draws what the canvas draws at its current zoom (zoom first to export a region)."
    />
  );
}

export const SvgAndPng: Story = {
  name: "SVG and PNG export",
  render: () => (
    <ExportFrame>
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <Case14 />
      </Suspense>
      <Genome />
    </ExportFrame>
  ),
};
