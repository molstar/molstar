/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `groups` behavior's body (_spec/README.md "Behaviors", Groups: lazy), loaded when a spec holds a section with a
// group (`groups.ts`, the door): a section's track groups folded as session state, the clade folds' model. A toggle is a hand's move (no step); collapse all, expand all and a reveal are jumps (back restores
// the folds with the window); `groups.keep` writes the folds shown into the groups' `collapsed` (one patch), and a
// patch or an undo that changes a group's `collapsed` resets the section's folds to the spec's.

import type { Op } from '@molstar/core/repr/patch';
import type { Node, Value } from '@molstar/core/repr/tree';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { EditOutcome, SpecModel, ViewModel } from '@molstar/seqx-viewer/models';
import type { ResolvedSection } from '@molstar/seqx-core/view';
import { GroupsCommandIds, type GroupsBehavior, type GroupsContext } from './groups.js';

const idOf = (address: string) => address.slice(address.lastIndexOf('/') + 1);
const sectionOf = (address: string) => address.slice(0, address.lastIndexOf('/'));

/** The commands, the reveal's opening and the reset of the folds on an edit of `collapsed`. */
export function groupsBody(spec: SpecModel, view: ViewModel, say?: (text: string) => void): GroupsBehavior {
  let disposed = false;
  const sections = (): readonly ResolvedSection[] => view.resolved.value?.sections ?? [];
  const section = (address: string | undefined) => sections().find((s) => s.address === address);
  const groupOf = (address: string | undefined) =>
    address === undefined ? undefined : section(sectionOf(address))?.trackGroups?.find((g) => g.address === address);
  /** The sections a bulk command acts on: the one named, else every section with groups. */
  const bulk = (ctx: GroupsContext) =>
    sections().filter((s) => s.trackGroups && (ctx.section === undefined || s.address === ctx.section));
  /** A command's `when`: the viewer live, then what it needs (true, or why not). */
  const when =
    (needs: (ctx: GroupsContext) => boolean, why: string) =>
    (ctx: GroupsContext): true | string =>
      disposed ? 'The viewer is retired' : needs(ctx) ? true : why;

  /** Every named section's folds set to all or none, as one jump (back restores them with the window). */
  const all = (ctx: GroupsContext, fold: boolean, done: string): CommandResult => {
    const o = view.jump(
      () => {
        let moved = false;
        for (const s of bulk(ctx)) {
          const now = view.foldedGroups(s.address);
          const want = fold ? s.trackGroups!.map((g) => idOf(g.address)) : [];
          if (want.length === now.length && want.every((id) => now.includes(id))) continue;
          moved = view.setFoldedGroups(s.address, want) || moved;
        }
        return moved;
      },
      true,
      'move',
    );
    return { message: o === 'moved' ? done : 'Already so' };
  };

  const said =
    (done: string) =>
    (o: EditOutcome): CommandResult =>
      o.ok
        ? { message: o.changed ? done : 'Already so', undoable: o.changed }
        : { message: o.issues.find((i) => i.severity === 'error')?.message ?? 'Refused' };

  const commands: CommandDefinition<GroupsContext>[] = [
    {
      id: GroupsCommandIds.toggle,
      title: 'Fold or open a group of tracks',
      category: 'Groups',
      when: when((ctx) => !!groupOf(ctx.address), 'No group there'),
      run: (_, ctx): CommandResult => {
        const g = groupOf(ctx.address)!;
        const s = sectionOf(g.address);
        const now = view.foldedGroups(s);
        const id = idOf(g.address);
        const folded = now.includes(id);
        // a hand's move: no step of the history (as a clade's toggle)
        view.setFoldedGroups(s, folded ? now.filter((x) => x !== id) : [...now, id]);
        // the announcer says the new state ("group Variants, collapsed, 3 tracks")
        say?.(`group ${g.label}, ${folded ? 'expanded' : 'collapsed'}, ${g.count} track${g.count === 1 ? '' : 's'}`);
        return { message: folded ? `${g.label} opened` : `${g.label} folded` };
      },
    },
    {
      id: GroupsCommandIds.collapseAll,
      title: 'Collapse all groups',
      category: 'Groups',
      when: when((ctx) => bulk(ctx).length > 0, 'No groups there'),
      run: (_, ctx) => all(ctx, true, 'Groups folded'),
    },
    {
      id: GroupsCommandIds.expandAll,
      title: 'Expand all groups',
      category: 'Groups',
      when: when((ctx) => bulk(ctx).length > 0, 'No groups there'),
      run: (_, ctx) => all(ctx, false, 'Groups opened'),
    },
    {
      id: GroupsCommandIds.keep,
      title: 'Keep these folds as the opening state',
      category: 'Groups',
      when: when((ctx) => !!section(ctx.section)?.trackGroups, 'No groups there'),
      run: (_, ctx) => {
        const s = section(ctx.section)!;
        const folded = view.foldedGroups(s.address);
        const ops: Op[] = s.trackGroups!.map((g) => {
          const at = g.address.split('/').filter(Boolean);
          return folded.includes(idOf(g.address))
            ? { op: 'set', at, key: 'collapsed', value: true }
            : { op: 'unset', at, key: 'collapsed' };
        });
        const o = spec.patch('Keep the folds', (_root, revision) => ({ base: revision, ops }), {
          ...(ctx.origin ? { origin: ctx.origin } : {}),
        });
        return o instanceof Promise ? o.then(said('Kept')) : said('Kept')(o);
      },
    },
  ];

  // a patch or an undo that changes a group's `collapsed` (or adds a group with `collapsed: true`) resets its section's
  // folds to the spec's; a group gone leaves the session's set, and one renamed (one gone, one new) keeps its fold
  const seen = new Map<string, Map<string, Value | undefined>>();
  const scan = (root: Node, tell: boolean) => {
    const now = new Map<string, Map<string, Value | undefined>>();
    const walk = (n: Node, address: string) => {
      for (const c of n.children ?? []) {
        if (c.id === undefined || c.kind === 'track' || c.kind === 'data') continue;
        if (c.kind === 'group') {
          let m = now.get(address);
          if (!m) now.set(address, (m = new Map()));
          m.set(c.id, c.params?.collapsed);
        } else walk(c, `${address}/${c.id}`);
      }
    };
    walk(root, '');
    const before = new Map(seen);
    seen.clear();
    for (const [k, v] of now) seen.set(k, v);
    if (!tell) return;
    for (const [s, own] of view.groupFolds.value) {
      const was = before.get(s) ?? new Map<string, Value | undefined>();
      const is = now.get(s) ?? new Map<string, Value | undefined>();
      let reset = false;
      const appeared: string[] = [];
      for (const [id, c] of is) {
        if (!was.has(id)) {
          appeared.push(id);
          if (c === true) reset = true;
        } else if (was.get(id) !== c) reset = true;
      }
      if (reset) {
        // the session's set already what the spec now opens with (`groups.keep`'s own patch): kept, nothing changes
        const opening = collapsedIds(root, s);
        if (own.length === opening.length && own.every((id) => opening.includes(id))) continue;
        view.setFoldedGroups(s, undefined);
        continue;
      }
      const gone = [...was.keys()].filter((id) => !is.has(id));
      if (!gone.length) continue;
      const kept = own.filter((id) => is.has(id));
      if (gone.length === 1 && appeared.length === 1 && own.includes(gone[0])) kept.push(appeared[0]);
      view.setFoldedGroups(s, kept);
    }
  };
  scan(spec.model.value.spec, false);
  const stop = spec.model.listen((model) => scan(model.spec, true));

  return {
    commands,
    open(tracks) {
      let moved = false;
      for (const t of tracks) {
        const s = section(sectionOf(sectionOf(t)));
        const g = s?.trackGroups?.find((x) => t.startsWith(`${x.address}/`));
        if (!g || !view.isFolded(g.address)) continue;
        const now = view.foldedGroups(s!.address);
        moved =
          view.setFoldedGroups(
            s!.address,
            now.filter((x) => x !== idOf(g.address)),
          ) || moved;
      }
      return moved;
    },
    dispose() {
      disposed = true;
      stop();
    },
  };
}

/** The ids of a section's groups the spec opens folded (`collapsed`), read from the spec itself. */
function collapsedIds(root: Node, section: string): string[] {
  let n: Node | undefined = root;
  for (const id of section.split('/').filter(Boolean)) n = n?.children?.find((c) => c.id === id);
  return (n?.children ?? []).flatMap((c) =>
    c.kind === 'group' && c.params?.collapsed === true && c.id !== undefined ? [c.id] : [],
  );
}
