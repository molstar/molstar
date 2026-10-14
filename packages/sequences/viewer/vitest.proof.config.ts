/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The node tier with the composer's proof (_test/proof/composer-frames.ts): `SEQ_PROOF_OUT=<folder>` writes every
// composer frame's hashes there, to compare a change that must draw as before.
import { defineConfig } from "vitest/config";
import { nodeTests } from "../vitest.shared.ts";

const base = nodeTests(import.meta.dirname, "all");
// Every frame recorded, under load too: a test timed out leaves its frames out and the folders differ (LC-45).
export default defineConfig({
  test: { ...base.test, setupFiles: ["_test/proof/composer-frames.ts"], testTimeout: 180_000, hookTimeout: 180_000 },
});
