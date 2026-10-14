/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a host's Storybook (.storybook/main.ts) takes from this package: where the stories and docs pages are, and a
// Vite plugin. The host also imports `@molstar/seq-storybook/tokens.css` in its preview.

import { fileURLToPath } from "node:url";

const src = fileURLToPath(new URL("./src", import.meta.url));

/** The globs of the package's stories and docs pages, for the host's `stories`. */
export const seqstarStories = [`${src}/_stories/docs/**/*.mdx`, `${src}/_stories/**/*.stories.@(ts|tsx)`];

/**
 * The Seq* stories compile the format's acceptance cases on the format's test fixtures, which read the Seq* JSON
 * fixtures (lib/seqstar/core/_test/fixtures) with `node:fs`: for those test modules only, `node:fs` is a small module serving
 * the bundled files (src/_stories/support/seqstar-fixture-fs.ts).
 */
export function seqstarFixtureFs() {
  const shim = `${src}/_stories/support/seqstar-fixture-fs.ts`;
  return {
    name: "seqstar-fixture-fs",
    enforce: "pre" as const,
    resolveId(id: string, importer?: string) {
      return id === "node:fs" && importer && /\/lib\/seqstar\/[^/]+\/_test\//.test(importer) ? shim : null;
    },
  };
}
