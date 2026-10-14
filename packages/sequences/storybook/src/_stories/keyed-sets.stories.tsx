/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useSeqViewer } from "@molstar/seq-viewer/react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { useMemo, useState } from "react";
import pf00042 from "../../../core/_test/fixtures/pf00042.json?raw";
import { type FamilyRow, keyedFamilyDoc, MODES, type Mode, modePatch } from "../../../viewer/_test/keyed-family.js";
import f from "./support/seqstar-format.module.css";
import { storyPatch } from "./support/seqstar-patch.js";
import { PrismPage, prism } from "./support/seqstar-prism.js";
import { ViewerCard } from "./support/seqstar-viewer.js";
import { Choices } from "./ui/index.js";

/**
 * Keyed sets: the members of a set are drawn in an order and kept by a
 * filter, both over key tables, `rows` annotations that give members fields. An alignment's key tables here are the
 * kernel providers `row_stats` (identity to a reference row, coverage, length) and `redundancy` (greedy clustering at
 * an identity threshold, in a given order: each row's `cluster` and whether it is its `representative`), computed at
 * load from the alignment and recomputed when it changes. A `rows` layer's `order` takes `"tree"`, `"data"` or a key
 * table's field (several, lexicographic: cluster, then identity), and `filter` keeps rows by a field's values or range.
 */
const meta = { title: "Seq*/Views and layouts/Alignments", parameters: { layout: "fullscreen" } } satisfies Meta;
export default meta;
type Story = StoryObj;

const FAMILY: FamilyRow[] = (JSON.parse(pf00042) as { rows: [string, number, string][] }).rows.map(
  ([label, , gapped]) => ({
    label,
    gapped,
  }),
);

/** The choices: each patches the `members` layer's `order` and `filter`, as the Stories agent would. */
function Modes() {
  const { context } = useSeqViewer();
  const [mode, setMode] = useState("data");
  const pick = (key: string) => {
    const m = MODES.find((x) => x.key === key) as Mode;
    setMode(m.key);
    void storyPatch(context, modePatch(m) as never, m.label);
  };
  return (
    <div className={prism.wide}>
      <Choices
        label="Order and filter"
        choices={MODES.map((m) => ({ id: m.key, label: m.label }))}
        active={mode}
        onPick={pick}
      />
    </div>
  );
}

export const OnAlignment: Story = {
  name: "Order and filter rows",
  render: () => {
    const d = useMemo(() => keyedFamilyDoc(FAMILY), []);
    return (
      <PrismPage
        title="Rows by identity and by cluster"
        caption={
          <>
            PF00042's 32 globins as the repository commits them. <code>row_stats</code> gives each row its identity to
            the reference row (HBA_HUMAN) and its coverage; <code>redundancy</code> takes the rows from the most to the
            least identical to the reference and puts each in the cluster of the first representative it is at least 90%
            (or 50%) identical to. All four are data columns. Pick an order or filter: they patch the layer's{" "}
            <code>order</code> and <code>filter</code> (key tables are read by field name). The family is diverse, so at
            90% only two pairs merge (30 representatives of 32); at 50% 22 rows represent the rest. "Representatives"
            keeps one row per cluster, the vertical map and the lanes following the kept rows.
          </>
        }
      >
        <div className={f.stack}>
          <ViewerCard doc={d} label="PF00042" options={{ maxHeight: 480 }}>
            <Modes />
          </ViewerCard>
        </div>
      </PrismPage>
    );
  },
};
