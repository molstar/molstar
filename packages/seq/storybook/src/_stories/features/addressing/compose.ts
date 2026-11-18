/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/core/repr/tree';
import type { SeqPorts } from '@molstar/seqx-core/context';
import type { ExampleView } from '@molstar/seqx-examples';
import type { LinkRun } from '@molstar/seqx-viewer/behaviors';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import type { RefEntry } from '@molstar/seqx-viewer/models';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';

// The addressing stories' composition, DOM-free: the story and its headless twin call these, so the twin composes the
// viewer and acts on it as the story does. Every action goes through the viewer's commands or its address model (a
// link through `runLink`), never around them.

export interface ComposeOptions {
  /** Where the canvas mounts: none in the headless twin. */
  readonly target?: HTMLElement;
  /** The delay of the view's delayed sources (the view's own pacing when absent). */
  readonly delayMs?: number;
  /** The ports: the twin's reads the fixtures from memory. */
  readonly ports?: SeqPorts;
}

/** A viewer of an example view, composed from the default spec with the view's registry. */
export function composeViewer(view: ExampleView, options: ComposeOptions = {}): SeqViewer {
  const base = DefaultSeqViewerSpec();
  const extra = view.registry?.(options.delayMs !== undefined ? { delayMs: options.delayMs } : {}) ?? [];
  return createSeqViewer(
    {
      ...base,
      ...(options.target ? {} : { ui: undefined }),
      registry: [...(base.registry ?? []), ...extra],
      ...(options.ports ? { ports: options.ports } : {}),
    },
    options.target,
  );
}

/** A bookmark as the list shows it: its label, what it names, and whether it resolves now. */
export interface BookmarkRow {
  readonly address: string;
  readonly label: string;
  readonly color?: string;
  readonly to: string;
  readonly status: RefEntry['status'];
}

export function bookmarkRows(viewer: SeqViewer): BookmarkRow[] {
  return viewer.addresses.refs.value.map((r) => ({
    address: r.address,
    label: r.label,
    ...(r.color !== undefined ? { color: r.color } : {}),
    to: r.to.join(' '),
    status: r.status,
  }));
}

/** Focuses a bookmark through the `address.focus` command: its message (a waiting one focuses when it resolves). */
export function focusBookmark(viewer: SeqViewer, address: string): string {
  return viewer.runCommand('address.focus', { address });
}

/**
 * Bookmarks what the view shows now (the region drawn, as an address: never another ref, so bookmarks never chain)
 * through `refs.add`: one patch, one step of the history. The id is the next free `view-<n>`; the result is the command's
 * message.
 */
export function bookmarkView(viewer: SeqViewer): string | Promise<string> {
  const to = viewer.addresses.shownAddress();
  if (to === undefined) return 'Nothing shown';
  const taken = new Set(viewer.addresses.refs.value.map((r) => r.id));
  let n = 1;
  while (taken.has(`view-${n}`)) n++;
  // through the history when the viewer has one: then it runs in turn
  return viewer.execute('refs.add', { ref: { id: `view-${n}`, to, label: `View ${n}` } });
}

/** Runs a link's actions: what a click on a link in the text does. */
export function followLink(viewer: SeqViewer, href: string): LinkRun | undefined {
  return viewer.runLink(href);
}

/** A link's hover: each address it names highlighted on the links' hover channel (EC-182); `leave` ends it. */
export function hoverLink(viewer: SeqViewer, href: string): void {
  const parsed = viewer.addresses.parseLink(href);
  for (const a of parsed?.link?.actions ?? []) viewer.runCommand('address.hover', { address: a.address });
}

export function leaveLink(viewer: SeqViewer): void {
  viewer.runCommand('address.unhover', {});
}

/** What the status line says of a run: its messages, then its issues' messages. */
export function describeRun(run: LinkRun | undefined): string {
  if (!run) return 'Not a link';
  const issues: readonly Issue[] = run.issues;
  return [...run.outcomes.map((o) => o.message), ...issues.map((i) => i.message)].filter(Boolean).join('; ');
}
