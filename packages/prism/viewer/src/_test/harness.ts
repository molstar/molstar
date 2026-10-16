/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { afterEach, expect } from 'vitest';
import type { ApiResult } from '../api/area.js';
import type { RunRecord } from '../api/run.js';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { createPrismViewerHeadless } from '../headless.js';
import type { PrismViewerSpec, ViewerExtension } from '../spec.js';
import type { PrismViewer } from '../viewer.js';

// A headless viewer for tests (viewer-spec.md §6): made per test, disposed after it, with the repository's example files.
// Tests never reach the network.

const dataDir = path.resolve(import.meta.dirname, '../../../../../data/examples');

/** An example file from `data/examples`, as a DOM `File` (named `as`, the format is read from the name). */
export function exampleFile(name: string, as = name): File {
  return new File([fs.readFileSync(path.join(dataDir, name))], as);
}

/** 1CBS (protein, retinoic acid ligand, waters): the usual fixture. */
export const cbs = () => exampleFile('1cbs_updated.cif', '1cbs.cif');
/** 1CRN (crambin, protein only). */
export const crambin = () => exampleFile('1crn.cif');

// Every viewer a test made, disposed after it (a test may make several).
const created: PrismViewer[] = [];
afterEach(async () => {
  const viewers = created.splice(0);
  for (const v of viewers.reverse()) await v.dispose();
});

/** Waits until `check` holds (polled), failing after `ms`: for the end of a fire-and-forget run, not a fixed sleep. */
export async function until(check: () => boolean, ms = 3000): Promise<void> {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms waiting for a condition`);
    await new Promise((r) => setTimeout(r, 5));
  }
}

/** The runs that finish from now on (live list): wait for one with `until(() => runs.some(…))`. */
export function recordRuns(v: PrismViewer): RunRecord[] {
  const runs: RunRecord[] = [];
  v.api.events.finished.subscribe((r) => runs.push(r));
  return runs;
}

/** The default list with these added; one with a default's id takes its place (a feature's tests pass the feature). */
export function withExtensions<E extends { readonly id: string }>(base: readonly E[], added: readonly E[]): E[] {
  const ids = new Set(added.map((e) => e.id));
  return [...base.filter((e) => !ids.has(e.id)), ...added];
}

/** A headless viewer with the default spec (and `spec` over it, `extensions` added), disposed after the test; `files` loaded first. */
export async function headlessViewer(
  o: { extensions?: ViewerExtension[]; files?: File[]; spec?: Partial<PrismViewerSpec> } = {},
): Promise<PrismViewer> {
  const base = DefaultPrismViewerSpec();
  const viewer = await createPrismViewerHeadless({
    spec: { ...base, ...o.spec, extensions: withExtensions(base.extensions ?? [], o.extensions ?? []) },
  });
  created.push(viewer);
  for (const f of o.files ?? []) expectOk(await viewer.api.load.files({ files: [f] }));
  return viewer;
}

/** Asserts a result is ok and narrows it. */
export function expectOk<R>(r: ApiResult<R>): asserts r is Extract<ApiResult<R>, { ok: true }> {
  if (!r.ok) expect.fail(`${r.code}: ${r.message}${r.hint ? ` (${r.hint})` : ''}`);
}

/** The representation types on the scene, sorted. */
export function representationTypes(v: PrismViewer): string[] {
  return v.plugin.managers.structure.hierarchy.current.structures
    .flatMap((s) =>
      s.components.flatMap((c) => c.representations.map((r) => (r.cell.transform.params as any).type.name as string)),
    )
    .sort();
}
