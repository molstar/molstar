/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { findNames, type SourceNameHit } from '@molstar/seqx-core/data';
import { rangesLocation } from '@molstar/seqx-core/query';
import type { SeqViewer } from './viewer.js';

// The locus box's reading (_spec/README.md "Commands", A locus): a host's box calls `goToLocus` with what the reader
// typed. The text is a locus first (`chr17:7,668,402-7,687,550`, `17:7668402`, a sequence's name or alias: the
// `(locus :text …)` form under the section, through the assembly's names rule); a text no sequence is named falls back
// to the sources' names (`findNames`: a GFF3's `Name`, `ID`, `Alias`, EC-228), its first hit a section shows focused
// (one on a sequence no section shows is said, not focused); an aborted search settles as `cancelled`.

/**
 * What a locus box's text went to: a locus, a name a source found (focused), a name found only on a sequence no section
 * shows (not focused, and why), the name search aborted (`cancelled`), or nothing (and why).
 */
export type LocusOutcome =
  | { readonly kind: 'locus' }
  | { readonly kind: 'name'; readonly hit: SourceNameHit }
  | { readonly kind: 'hidden'; readonly hit: SourceNameHit; readonly message: string }
  | { readonly kind: 'cancelled' }
  | { readonly kind: 'none'; readonly message: string };

/** Quotes a text for a form's string (`"…"`, its quotes and backslashes escaped). */
function quoted(text: string): string {
  return `"${text.replace(/\\/g, '\\\\').replace(/"/g, '\\"')}"`;
}

/** The first section of assemblies of the view (its `of` one or a list of assemblies), when no section is named. */
function assemblySection(viewer: SeqViewer): string | undefined {
  const model = viewer.spec.model.value;
  const root = viewer.spec.shown.value;
  for (const s of viewer.view.resolved.value?.sections ?? []) {
    let n = root;
    for (const id of s.address.split('/').filter(Boolean)) n = n.children?.find((c) => c.id === id) ?? n;
    const of = n.params?.of;
    const list = (Array.isArray(of) ? of : [of]).filter((x): x is string => typeof x === 'string');
    if (list.length && list.every((a) => model.assemblies.has(a))) return s.address;
  }
  return undefined;
}

/**
 * Goes to what a locus box's text names (_spec/README.md "Commands", A locus): a locus under `section` (the first
 * section of assemblies when absent), else the first item a source finds by that name, its whole extent focused (a
 * jump, a step of the history). `signal` aborts the name search.
 */
export async function goToLocus(
  viewer: SeqViewer,
  text: string,
  options: { readonly section?: string; readonly signal?: AbortSignal } = {},
): Promise<LocusOutcome> {
  const t = text.trim();
  if (!t) return { kind: 'none', message: 'Nothing typed' };
  const section = options.section ?? assemblySection(viewer);
  if (section) {
    const { resolution, location } = viewer.addresses.locate(`${section}/(locus :text ${quoted(t)})`);
    if (location.parts.length) {
      viewer.addresses.focusLocation(location);
      return { kind: 'locus' };
    }
    const unknown = resolution.issues.some((i) => i.code === 'seq.query.locus.name');
    if (!unknown) return { kind: 'none', message: resolution.issues[0]?.message ?? 'Not a locus' };
  }
  // no sequence of that name: a source's item by name (its system named by the source, read through the assembly's
  // names rule onto its sequence)
  const signal = options.signal ?? new AbortController().signal;
  let hits: readonly SourceNameHit[];
  try {
    hits = await findNames(viewer.store.store.named(), t, { max: 20 }, signal);
  } catch (e) {
    // a box typed on (the host aborted the search): settled as cancelled, never an error to the caller
    if (signal.aborted) return { kind: 'cancelled' };
    throw e;
  }
  if (signal.aborted) return { kind: 'cancelled' };
  const model = viewer.spec.model.value;
  // the first hit a section shows is focused; a hit on a sequence no section shows is said, never focused unseen
  let hidden: { hit: SourceNameHit; system: string } | undefined;
  for (const hit of hits) {
    const system = model.systems.has(hit.system) ? hit.system : sequenceNamed(viewer, hit.system);
    if (!system) continue;
    if (!shownSystem(viewer, system)) {
      hidden ??= { hit, system };
      continue;
    }
    viewer.addresses.focusLocation(rangesLocation(system, [hit.start, hit.end]));
    return { kind: 'name', hit };
  }
  if (hidden) {
    const on = model.systems.get(hidden.system)?.name ?? hidden.hit.system;
    return { kind: 'hidden', hit: hidden.hit, message: `${t} is on ${on}, which no section shows` };
  }
  return { kind: 'none', message: `Nothing is named ${t}` };
}

/** Whether a section shows a system: its axis's, or a block's (a merged copy's too). */
function shownSystem(viewer: SeqViewer, system: string): boolean {
  for (const s of viewer.view.resolved.value?.sections ?? []) {
    if (s.axis.system === system) return true;
    if (s.blocks?.some((b) => b.axis.system === system || b.systems.includes(system))) return true;
  }
  return false;
}

/** The derived sequence a source's system name is, through each assembly's names rule; undefined when none. */
function sequenceNamed(viewer: SeqViewer, name: string): string | undefined {
  const model = viewer.spec.model.value;
  for (const a of model.assemblies.keys()) {
    const s = model.assemblySequence?.(a, name);
    if (s && s !== 'pending') return s.address;
  }
  return undefined;
}
