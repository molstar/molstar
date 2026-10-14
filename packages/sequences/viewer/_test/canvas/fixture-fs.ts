/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `node:fs` for the Seq* test fixtures in the browser (vitest.browser.config.ts sends their `node:fs` import here):
// the format's load-fixtures.ts and the model's alignment-fixtures.ts read the Seq* JSON fixtures (lib/seqstar/core/_test/fixtures)
// with `readFileSync(new URL(…, import.meta.url))`; here those files come bundled, by name. Nothing else is served.
// The same shim as Storybook's (lib/seqstar/storybook/src/_stories/support/seqstar-fixture-fs.ts).

import brs from "../../../core/_test/fixtures/1brs.json?raw";
import iep from "../../../core/_test/fixtures/1iep.json?raw";
import pae from "../../../core/_test/fixtures/af-barnase-barstar-pae.json?raw";
import p53 from "../../../core/_test/fixtures/p53.json?raw";
import pf00042 from "../../../core/_test/fixtures/pf00042.json?raw";

const FILES: Readonly<Record<string, string>> = {
  "1iep.json": iep,
  "1brs.json": brs,
  "af-barnase-barstar-pae.json": pae,
  "p53.json": p53,
  "pf00042.json": pf00042,
};

export function readFileSync(path: string | URL): string {
  const name = String(path instanceof URL ? path.pathname : path)
    .split("/")
    .pop() as string;
  const text = FILES[name];
  if (text === undefined) throw new Error(`no bundled fixture ${name}`);
  return text;
}
