/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from "@storybook/react-vite";
import { useMemo } from "react";
import { circularGenomeDoc } from "./support/seqstar-docs-circular.js";
import { PrismPage } from "./support/seqstar-prism.js";
import { ViewerCard, type ViewerOptions } from "./support/seqstar-viewer.js";

/**
 * Circular maps (CG2, CG3): a circular genome through the builder's preset, rings with their linear detail and its
 * window (the plasmid map is case 16, in the Viewer cases).
 */
const meta = { title: "Seq*/Views and layouts/Axes and strands", parameters: { layout: "fullscreen" } } satisfies Meta;
export default meta;
type Story = StoryObj;

/** Clean chrome (`UI.Chrome`): no part of the default UI's controls is needed here. */
const CLEAN: ViewerOptions = { chrome: "clean" };

export const CircularGenome: Story = {
  name: "A circular genome",
  render: () => {
    const doc = useMemo(circularGenomeDoc, []);
    return (
      <PrismPage
        title="A circular genome: strands on rings, a values ring, its linear detail"
        caption={
          <>
            <code>b.preset.circularGenome(view, chr, {"{ features, values }"})</code>: 160 synthetic genes on 120,000
            nt, the forward strand&apos;s on the outer ring and the reverse strand&apos;s inside it (each with
            arrowheads where a gene is long enough to carry one), a GC skew as bars on the inner ring, and the linear
            detail following, its view a window on the ring: wheel-zoom the detail, then drag the window or its grips.
          </>
        }
      >
        <ViewerCard doc={doc} label="chr, ring and detail" options={CLEAN} />
      </PrismPage>
    );
  },
};
