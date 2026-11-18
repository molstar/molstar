/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `editing` behavior (_spec/README.md "Behaviors"): edits of the spec as commands, each one patch through
// `SpecModel.patch` (one step of the host's history), so a menu, a key, an agent and a test edit alike: hide or show a
// layer or a whole track, move a track among its siblings, set a layer's color (its theme), remove a ref; story steps
// entered and left (navigation); and a gesture's preview, a model method and never a command (it runs per pointer
// move), that ends in one patch (`endPreview`) or nothing (`cancelPreview`, Esc).

import type { Op, Patch } from '@molstar/core/repr/patch';
import { type Node, resolveAddress, type Value } from '@molstar/core/repr/tree';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { EditOutcome, SpecModel } from '@molstar/seqx-viewer/models';
import type { EditContext, EditOrigin } from '@molstar/seqx-viewer/ports';

/** What an editing command is given: the node it edits, and what each needs. */
export interface EditingContext {
  /** The node's address: a layer (`layer.*`), a track (`track.*`), a ref (`refs.remove`). */
  readonly address?: string;
  /** `track.move`: the sibling it goes after; null: first. */
  readonly after?: string | null;
  /** `layer.color`: the theme form or color the layer is drawn with. */
  readonly color?: Value;
  /** `step.enter`: the story step's id. */
  readonly step?: string;
  /** Who edits (`person` when absent; an agent's turn joins one step). */
  readonly origin?: EditOrigin;
  /** The host's run the edit is part of (a host's step that runs this command). */
  readonly within?: EditContext;
}

/** The editing commands' ids. */
export const EditingCommandIds = {
  hideLayer: 'layer.hide',
  showLayer: 'layer.show',
  colorLayer: 'layer.color',
  hideTrack: 'track.hide',
  showTrack: 'track.show',
  moveTrack: 'track.move',
  removeRef: 'refs.remove',
  enterStep: 'step.enter',
  leaveStep: 'step.leave',
  cancelPreview: 'preview.cancel',
} as const;

export interface EditingBehavior {
  readonly commands: readonly CommandDefinition<EditingContext>[];
  /** A gesture's preview: a patch over the root (a drag resizing a track), replacing the last; false when refused. */
  preview(patch: Patch): boolean;
  /** Ends the gesture in its one patch: one step of the history, or nothing when refused or there is none. */
  endPreview(label: string, origin?: EditOrigin): EditOutcome | Promise<EditOutcome> | undefined;
  /** Ends the gesture in nothing (Esc). */
  cancelPreview(): void;
  dispose(): void;
}

const ids = (address: string) => address.split('/').filter(Boolean);

/** The node at an address of ids, when it is of the kind. */
function nodeOf(root: Node, address: string | undefined, kind: string): Node | undefined {
  if (address === undefined) return undefined;
  const n = resolveAddress(root, ids(address));
  return n?.kind === kind ? n : undefined;
}

export function editingBehavior(spec: SpecModel): EditingBehavior {
  let disposed = false;
  const said =
    (done: string) =>
    (o: EditOutcome): CommandResult =>
      o.ok
        ? { message: o.changed ? done : 'Already so', undoable: o.changed }
        : { message: o.issues.find((i) => i.severity === 'error')?.message ?? 'Refused' };
  /** One patch made when the edit runs, against the revision it runs on: commands issued together each apply. */
  const edit = (label: string, ops: (root: Node) => Op[], ctx: EditingContext, done: string) => {
    const o = spec.patch(label, (root, revision) => ({ base: revision, ops: ops(root) }), {
      ...(ctx.origin ? { origin: ctx.origin } : {}),
      ...(ctx.within ? { within: ctx.within } : {}),
    });
    return o instanceof Promise ? o.then(said(done)) : said(done)(o);
  };
  const needs =
    (kind: string, what: string) =>
    (ctx: EditingContext): true | string =>
      disposed ? 'The viewer is retired' : nodeOf(spec.root.value, ctx.address, kind) ? true : `No ${what} there`;

  const layersOf = (root: Node, track: string) =>
    (nodeOf(root, track, 'track')?.children ?? []).flatMap((c) => (c.id === undefined ? [] : [c.id]));

  const commands: CommandDefinition<EditingContext>[] = [
    {
      id: EditingCommandIds.hideLayer,
      title: 'Hide a layer',
      category: 'Edit',
      when: needs('layer', 'layer'),
      run: (_, ctx) =>
        edit('Hide a layer', () => [{ op: 'set', at: ids(ctx.address!), key: 'hidden', value: true }], ctx, 'Hidden'),
    },
    {
      id: EditingCommandIds.showLayer,
      title: 'Show a layer',
      category: 'Edit',
      when: needs('layer', 'layer'),
      run: (_, ctx) =>
        edit('Show a layer', () => [{ op: 'unset', at: ids(ctx.address!), key: 'hidden' }], ctx, 'Shown'),
    },
    {
      id: EditingCommandIds.colorLayer,
      title: 'Color a layer',
      category: 'Edit',
      when: (ctx) => {
        const can = needs('layer', 'layer')(ctx);
        return can !== true ? can : ctx.color === undefined ? 'No color' : true;
      },
      run: (_, ctx) =>
        edit(
          'Color a layer',
          () => [{ op: 'set', at: ids(ctx.address!), key: 'color', value: ctx.color! }],
          ctx,
          'Colored',
        ),
    },
    {
      id: EditingCommandIds.hideTrack,
      title: 'Hide a track',
      category: 'Edit',
      when: needs('track', 'track'),
      run: (_, ctx) =>
        edit(
          'Hide a track',
          (root) =>
            layersOf(root, ctx.address!).map((l) => ({
              op: 'set',
              at: [...ids(ctx.address!), l],
              key: 'hidden',
              value: true,
            })),
          ctx,
          'Hidden',
        ),
    },
    {
      id: EditingCommandIds.showTrack,
      title: 'Show a track',
      category: 'Edit',
      when: needs('track', 'track'),
      run: (_, ctx) =>
        edit(
          'Show a track',
          (root) =>
            layersOf(root, ctx.address!).map((l) => ({ op: 'unset', at: [...ids(ctx.address!), l], key: 'hidden' })),
          ctx,
          'Shown',
        ),
    },
    {
      id: EditingCommandIds.moveTrack,
      title: 'Move a track',
      category: 'Edit',
      when: (ctx) => {
        const can = needs('track', 'track')(ctx);
        return can !== true ? can : ctx.after === undefined ? 'No place to move it to' : true;
      },
      run: (_, ctx) => {
        const at = ids(ctx.address!);
        return edit(
          'Move a track',
          () => [{ op: 'move', at, parent: at.slice(0, -1), place: { after: ctx.after ?? null } }],
          ctx,
          'Moved',
        );
      },
    },
    {
      id: EditingCommandIds.removeRef,
      title: 'Remove a bookmark',
      category: 'Edit',
      when: needs('ref', 'ref'),
      run: (_, ctx) => edit('Remove a bookmark', () => [{ op: 'remove', at: ids(ctx.address!) }], ctx, 'Removed'),
    },
    {
      id: EditingCommandIds.enterStep,
      title: 'Show a story step',
      category: 'Story',
      when: (ctx) =>
        disposed
          ? 'The viewer is retired'
          : ctx.step !== undefined && nodeOf(spec.root.value, `/snapshots/${ctx.step}`, 'snapshot')
            ? true
            : 'No step to show',
      run: (_, ctx): CommandResult => {
        const issue = spec.enterStep(ctx.step!);
        if (issue) return { message: issue.message };
        return { message: spec.step.value?.applies ? 'Shown' : 'The step no longer applies' };
      },
    },
    {
      id: EditingCommandIds.leaveStep,
      title: 'Leave the story step',
      category: 'Story',
      when: () => (disposed ? 'The viewer is retired' : spec.step.value ? true : 'No step shown'),
      run: (): CommandResult => {
        spec.leaveStep();
        return { message: 'Left' };
      },
    },
    {
      id: EditingCommandIds.cancelPreview,
      title: 'Cancel the gesture',
      category: 'Edit',
      shortcut: { key: 'Escape' },
      when: () => (disposed ? 'The viewer is retired' : spec.preview.value ? true : 'Nothing to cancel'),
      run: (): CommandResult => {
        spec.cancelPreview();
        return { message: 'Cancelled' };
      },
    },
  ];

  return {
    commands,
    preview: (patch) => !disposed && spec.setPreview(patch).ok,
    endPreview: (label, origin) => {
      if (disposed || !spec.preview.value) return undefined;
      return spec.edit(label, (e) => e.commitPreview(), origin ? { origin } : {});
    },
    cancelPreview: () => spec.cancelPreview(),
    dispose() {
      disposed = true;
      spec.cancelPreview();
    },
  };
}
