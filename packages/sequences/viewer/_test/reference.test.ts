/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's generated reference: written to the spec's reference folder, so a change to a command, a behavior, a
// config item or a built-in shows there.

import { describe, expect, it } from "vitest";
import { DefaultSeqViewerSpec, SeqCommands, SeqViewerConfig } from "../src/index.js";

import { printViewerReference } from "../src/reference.js";

describe("the viewer's generated reference", () => {
  it("matches lib/seqstar/_spec/reference/viewer.md", async () => {
    await expect(printViewerReference()).toMatchFileSnapshot("../../_spec/reference/viewer.md");
  });

  it("lists every built-in command with its description, every behavior and every config item", () => {
    const text = printViewerReference();
    for (const c of Object.values(SeqCommands)) {
      expect(c.description).not.toBe("");
      expect(text).toContain(`### \`${c.id}\``);
    }
    for (const { behavior } of DefaultSeqViewerSpec().behaviors) expect(text).toContain(`| \`${behavior.name}\` |`);
    for (const group of Object.values(SeqViewerConfig))
      for (const item of Object.values(group)) {
        expect(item.description).not.toBe("");
        expect(text).toContain(`\`${item.key}\``);
      }
  });
});
