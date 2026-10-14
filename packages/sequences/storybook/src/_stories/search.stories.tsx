/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from "@storybook/react-vite";
import { LocationSetsStory } from "./support/seqstar-location-sets.js";
import { SearchStory } from "./support/seqstar-search.js";

/**
 * Search is not in the format or the viewer: `@molstar/seq-core/algorithms/search` ships matchers,
 * an application composes its own search from them. This one has a box beside a mounted viewer, runs exact, degenerate
 * (IUPAC for DNA, a protein grouping the story defines) and approximate (Myers' bit-parallel infix search, up to 3 edits)
 * matchers over the p53 case's protein and the gene case's DNA and its translation, lights every hit through the viewer's
 * highlight channels, accents the current hit and scrolls the viewer to it, and steps with Next and Previous or Enter and
 * Shift+Enter. What a hit means (here: a range of the sequence's coordinate system) is the story's code.
 */
const meta = {
  title: "Seq*/Query and search/Search",
  parameters: { layout: "fullscreen", docs: { source: { type: "code" } } },
} satisfies Meta;
export default meta;
type Story = StoryObj<typeof meta>;

export const AppComposedSearch: Story = {
  name: "App-composed search",
  render: () => <SearchStory />,
};

/**
 * Location sets: a section's `highlights`, each entry a highlight channel
 * of the viewer's (`band` in the main pass, `highlight` and `accent` on the overlay), a SeqQL query node among them; and
 * bookmarks, an inline `ranges` annotation the `bookmark.add`, `bookmark.remove` and `bookmark.step` actions edit and step.
 */
export const LocationSets: Story = {
  name: "Location sets and bookmarks",
  render: () => <LocationSetsStory />,
};
