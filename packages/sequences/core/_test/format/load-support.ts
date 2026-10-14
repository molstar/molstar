/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Helpers of the load tests (load.test.ts).

import type { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { Loaders } from "../../src/format/index.js";

/** The document with a view, which every document needs. */
export const done = (b: ReturnType<typeof createSeqViewSpec>) => {
  b.view({ ref: "end", label: "V" }).track({ label: "L" }).letters();
  return b.getDocument();
};

/** A loader that serves `files` (by uri or host file) and fails everything else; `fail` names those that fail too. */
export function loaders(files: Record<string, string>, fail: (name: string) => boolean = () => false) {
  const requests: string[] = [];
  const l: Loaders & { requests: string[] } = {
    requests,
    async load(request) {
      const name = request.uri ?? request.hostFile ?? "";
      requests.push(name);
      if (fail(name) || files[name] === undefined) throw new Error(`503 ${name}`);
      return new TextEncoder().encode(files[name]);
    },
  };
  return l;
}
