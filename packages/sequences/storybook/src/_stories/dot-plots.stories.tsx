/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from "@storybook/react-vite";
import { useMemo } from "react";
import { calmodulinDotPlot } from "./support/seqstar-dot-plots.js";
import { PrismPage } from "./support/seqstar-prism.js";
import { ViewerCard } from "./support/seqstar-viewer.js";

/**
 * Dot plots: positions against positions on one section's axis, x across as every track's and the same window down
 * the track, so zooming and panning the tracks moves both. A grid's values are cells; pairs are segments.
 */
const meta = { title: "Seq*/Views and layouts/Dot plots", parameters: { layout: "fullscreen" } } satisfies Meta;
export default meta;
type Story = StoryObj;

export const SelfComparison: Story = {
  name: "A sequence against itself",
  render: function Render() {
    const doc = useMemo(calmodulinDotPlot, []);
    return (
      <PrismPage
        title="Calmodulin against itself: a duplication as parallel diagonals"
        caption={
          <>
            Calmodulin's two lobes are a duplication, each of two EF hands. The plot is the residues across against the
            same residues down: the identity over a window of 9, a third or more, as blue cells (the diagonal is the
            sequence against itself), and the 4-letter words found twice as red segments, off the diagonal by a lobe
            (about 72) and by an EF hand (about 36). Hover a cell: both residues light in the tracks above, with a
            crosshair through the cell. Drag on the plot to select a rectangle, its two ranges. Zoom or pan the tracks:
            the plot shows the same window both ways (Case 8 shows AlphaFold DB's predicted aligned error the same way).
          </>
        }
      >
        <ViewerCard doc={doc} label="Calmodulin (P0DP23) against itself" />
      </PrismPage>
    );
  },
};
