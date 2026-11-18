/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { printSexpr } from '@molstar/core/repr/codecs/sexpr';
import type { Patch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import type { SeqPorts } from '@molstar/seqx-core/context';
import type { ExampleView } from '@molstar/seqx-examples';
import type { EditingBehavior } from '@molstar/seqx-viewer/behaviors';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { HistoryParticipants, hostHistory } from '@molstar/seqx-viewer/ports';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { HistoryModel } from '@molstar/ui/history/history-model';

// The editing story's composition, DOM-free: the story and its headless twin call these. A stand-in host owns one undo
// stack (`HistoryModel` over `HistoryParticipants`) with two participants: the viewer's spec (through `hostHistory`)
// and a stand-in structure state the host's own command changes. Every edit goes through the viewer's commands or
// its editing behavior's preview, never around them.

export const Section = '/views/main/protein';
export const Tracks = ['letters', 'domains', 'conservation'] as const;
export const track = (id: string) => `${Section}/${id}`;

export interface EditingSession {
  readonly viewer: SeqViewer;
  /** The host's stack: Undo and Redo for both participants. */
  readonly stack: HistoryModel<Readonly<Record<string, unknown>>>;
  /** The stand-in structure's color: the host's own state, a participant of the same stack. */
  readonly structure: State<string>;
  dispose(): void;
}

export interface ComposeOptions {
  readonly target?: HTMLElement;
  readonly ports?: SeqPorts;
  /** The delay of the view's delayed sources (the view's own pacing when absent). */
  readonly delayMs?: number;
}

export function composeEditing(view: ExampleView, options: ComposeOptions = {}): EditingSession {
  const participants = new HistoryParticipants();
  const stack = new HistoryModel<Readonly<Record<string, unknown>>>(participants);
  const structure = new State('grey');
  const removeStructure = participants.add('structure', {
    checkpoint: () => structure.value,
    restore: async (c: string) => structure.set(c),
  });
  const base = DefaultSeqViewerSpec();
  const viewer = createSeqViewer(
    {
      ...base,
      ...(options.target ? {} : { ui: undefined }),
      registry: [
        ...(base.registry ?? []),
        ...(view.registry?.(options.delayMs !== undefined ? { delayMs: options.delayMs } : {}) ?? []),
      ],
      history: hostHistory(stack, participants),
      ...(options.ports ? { ports: options.ports } : {}),
    },
    options.target,
  );
  return {
    viewer,
    stack,
    structure,
    dispose() {
      viewer.retire();
      removeStructure();
      stack.dispose();
    },
  };
}

const editing = (s: EditingSession) => s.viewer.behaviors.get('editing') as EditingBehavior;

/** Whether a track is hidden (its own `hidden`). */
export function trackHidden(s: EditingSession, id: string): boolean {
  const t = s.viewer.spec.root.value.children
    ?.find((c) => c.id === 'views')
    ?.children?.[0]?.children?.find((c) => c.id === 'protein')
    ?.children?.find((c) => c.id === id);
  return t?.params?.hidden === true;
}

/** Hides or shows a whole track: one patch, one step. */
export function toggleTrack(s: EditingSession, id: string): Promise<string> {
  const command = trackHidden(s, id) ? 'track.show' : 'track.hide';
  return Promise.resolve(s.viewer.execute(command, { address: track(id) }));
}

/** Moves a track first among its siblings: one patch. */
export function moveFirst(s: EditingSession, id: string): Promise<string> {
  return Promise.resolve(s.viewer.execute('track.move', { address: track(id), after: null }));
}

/** Sets the domains' theme color: one patch. */
export function setDomainsColor(s: EditingSession, color: string): Promise<string> {
  return Promise.resolve(s.viewer.execute('layer.color', { address: `${track('domains')}/ranges`, color }));
}

/** Bookmarks the region shown: one patch. */
export function bookmarkView(s: EditingSession): Promise<string> {
  const to = s.viewer.addresses.shownAddress();
  if (to === undefined) return Promise.resolve('Nothing shown');
  const taken = new Set(s.viewer.addresses.refs.value.map((r) => r.id));
  let n = 1;
  while (taken.has(`view-${n}`)) n++;
  return Promise.resolve(s.viewer.execute('refs.add', { ref: { id: `view-${n}`, to, label: `View ${n}` } }));
}

/** The host's own step: the stand-in structure recolored. */
export async function recolorStructure(s: EditingSession, color: string): Promise<void> {
  await s.stack.run('Recolor the structure', { kind: 'person' }, async () => s.structure.set(color));
}

const heightPatch = (id: string, height: number): Patch => ({
  ops: [{ op: 'set', at: ['views', 'main', 'protein', id], key: 'height', value: height }],
});

/** A drag resizing a track: a preview each move, replacing the last. */
export function dragHeight(s: EditingSession, id: string, height: number): boolean {
  return editing(s).preview(heightPatch(id, height));
}

/** The drag's end: its one patch, one step. */
export async function endDrag(s: EditingSession): Promise<void> {
  await editing(s).endPreview('Resize a track');
}

/** Esc: the drag ends in nothing. */
export function cancelDrag(s: EditingSession): string {
  return s.viewer.runCommand('preview.cancel', {});
}

/** The stack's steps, newest first, and what Redo would put back. */
export function historyRows(s: EditingSession): { readonly done: string[]; readonly redo: string[] } {
  return {
    done: s.stack.steps.value.map((x) => x.label),
    redo: s.stack.redoable.value.map((x) => x.label),
  };
}

/** The spec shown (the preview included), without its data: the views and the refs, printed. */
export function printedSpec(root: Node): string {
  const kept = (root.children ?? []).filter((c) => c.id !== 'data');
  return printSexpr({ ...root, children: kept });
}
