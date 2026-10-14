/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `.svsx` story's data: the committed fixture opened with the
// reader and compiled over `archiveLoaders` (no network), and the bundle step run over a fixture fetch (offline).

import { bundleUris, writeSvsx } from "@molstar/seq-core/builder/svsx";
import { archiveLoaders, compile, readSvsx, type Svsx, toJSON } from "@molstar/seq-core/format";
import fixtureUrl from "../../../../core/_test/format/fixtures/kinase.svsx?url";
import { absoluteDoc, fixtureFetch, SVSX_FILES } from "../../../../core/_test/format/svsx-fixture.ts";
import { type CompiledCase, reg } from "./seqstar-cases";

const refuse = {
  async load(): Promise<Uint8Array> {
    throw new Error("the archive has no such file, and nothing else is asked");
  },
};

export interface Opened {
  readonly svsx: Svsx;
  readonly size: number;
  readonly compiled: CompiledCase | null;
}

/** The fixture archive, read and compiled over its own files: a loader that fails everything else is never asked. */
export async function openFixture(): Promise<Opened> {
  const bytes = new Uint8Array(await (await fetch(fixtureUrl)).arrayBuffer());
  const svsx = await readSvsx(bytes, { registries: reg });
  if (!svsx.doc) return { svsx, size: bytes.length, compiled: null };
  const { set, issues } = await compile(svsx.doc, reg, archiveLoaders(svsx.files, refuse, svsx.base));
  const views = (svsx.doc.root.children ?? [])
    .filter((c) => c.kind === "view")
    .map((c) => ({ ref: c.ref as string, label: String(c.params?.label ?? c.ref) }));
  return {
    svsx,
    size: bytes.length,
    compiled: { doc: svsx.doc, set, validation: svsx.issues, compile: issues, views },
  };
}

export interface Bundle {
  readonly before: string;
  readonly after: string;
  readonly requests: readonly string[];
  readonly sources: Readonly<Record<string, string>>;
  readonly read: Svsx;
  readonly size: number;
}

/** A document with absolute URIs bundled over a fixture fetch, written, and read back. */
export async function makeBundle(): Promise<Bundle> {
  const doc = absoluteDoc();
  const f = fixtureFetch();
  const bundled = await bundleUris(doc, f.fetch);
  const bytes = await writeSvsx(bundled.doc, bundled.files);
  return {
    before: toJSON(doc, { pretty: true }),
    after: toJSON(bundled.doc, { pretty: true }),
    requests: f.requests,
    sources: SVSX_FILES,
    read: await readSvsx(bytes, { registries: reg }),
    size: bytes.length,
  };
}
