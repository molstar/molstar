/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `groups` behavior (_spec/README.md "Behaviors", Groups): a section's track groups folded as session state, the
// clade folds' model. This is its door: the commands listed, their bodies (`groups-body.ts`) loaded the first time a
// spec holds a section with a group, as render's groups chunk (core's `withViewKind`: the model waits for both). Why
// lazy: groups are rare; a viewer that never shows one pays for little of them. `layer.lanes` is here whole: it acts
// on a ranges layer, grouped or not.

import type { Op } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { EditOutcome, SpecModel, ViewModel } from '@molstar/seqx-viewer/models';
import type { EditOrigin } from '@molstar/seqx-viewer/ports';
import type { BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { withViewKind } from '@molstar/seqx-core/resolve';

/** What a groups command is given: a group's address (`toggle`), a section's (`collapseAll`, `expandAll`, `keep`). */
export interface GroupsContext {
  readonly address?: string;
  readonly section?: string;
  readonly origin?: EditOrigin;
  /** `layer.lanes`: a ranges layer's lanes, `packed` or `one` (its `address` the layer's). */
  readonly lanes?: 'packed' | 'one';
}

export const GroupsCommandIds = {
  toggle: 'groups.toggle',
  collapseAll: 'groups.collapseAll',
  expandAll: 'groups.expandAll',
  keep: 'groups.keep',
  lanes: 'layer.lanes',
} as const;

export interface GroupsBehavior {
  readonly commands: readonly CommandDefinition<GroupsContext>[];
  /** The viewer composed: the announcer's notice for a fold. */
  composed?(viewer: SeqViewer): void;
  /** Opens the folded groups holding these tracks (a reveal landing on an item only they draw): no step of its own. */
  open(tracks: readonly string[]): boolean;
  dispose(): void;
}

let Body: typeof import('./groups-body.js') | undefined;
let loading: Promise<unknown> | undefined;
/** Loads the commands' bodies. */
const loadBody = () =>
  (loading ??= import('./groups-body.js').then(
    (m) => {
      Body = m;
    },
    () => {
      loading = undefined;
    },
  ));

/** Whether a section holds a group of tracks. */
function grouped(n: Node): boolean {
  const c = n.children;
  if (c) for (let i = 0; i < c.length; i++) if (c[i].kind === 'group') return true;
  return false;
}

withViewKind('section', loadBody, grouped);

const Commands: readonly (readonly [string, string, string])[] = [
  [GroupsCommandIds.toggle, 'Fold or open a group of tracks', 'No group there'],
  [GroupsCommandIds.collapseAll, 'Collapse all groups', 'No groups there'],
  [GroupsCommandIds.expandAll, 'Expand all groups', 'No groups there'],
  [GroupsCommandIds.keep, 'Keep these folds as the opening state', 'No groups there'],
];

export function groupsBehavior(spec: SpecModel, view: ViewModel): GroupsBehavior {
  let disposed = false;
  let body: GroupsBehavior | undefined;
  let viewer: SeqViewer | undefined;
  const say = (text: string) => viewer?.notify(text);
  // the body made once its chunk landed (the model that needed it is the next one): it reads the folds from then on
  const live = () => (body ??= !disposed && Body ? Body.groupsBody(spec, view, say) : undefined);
  live();
  const stop = spec.model.listen(() => void live());
  const said =
    (done: string) =>
    (o: EditOutcome): CommandResult =>
      o.ok
        ? { message: o.changed ? done : 'Already so', undoable: o.changed }
        : { message: o.issues.find((i) => i.severity === 'error')?.message ?? 'Refused' };
  const of = (id: string) => live()?.commands.find((c) => c.id === id);

  const commands: CommandDefinition<GroupsContext>[] = Commands.map(
    ([id, title, none]): CommandDefinition<GroupsContext> => ({
      id,
      title,
      category: 'Groups',
      when: (ctx: GroupsContext) => (disposed ? 'The viewer is retired' : (of(id)?.when?.(ctx) ?? none)),
      run: (args, ctx, signal) => of(id)!.run(args, ctx, signal),
    }),
  );
  commands.push({
    // a ranges layer in one lane (labels off) or packed: the figure's look, one patch (the layer's `lanes`)
    id: GroupsCommandIds.lanes,
    title: 'Lay a ranges layer in one lane or packed',
    category: 'Groups',
    when: (ctx) => {
      if (disposed) return 'The viewer is retired';
      if (ctx.address === undefined || (ctx.lanes !== 'one' && ctx.lanes !== 'packed')) return 'No ranges layer there';
      let n: Node | undefined = spec.root.value;
      for (const id of ctx.address.split('/').filter(Boolean)) n = n?.children?.find((c) => c.id === id);
      return n?.kind === 'layer' && n.params?.type === 'ranges' ? true : 'No ranges layer there';
    },
    run: (_, ctx) => {
      const at = ctx.address!.split('/').filter(Boolean);
      const ops: Op[] = [
        ctx.lanes === 'one' ? { op: 'set', at, key: 'lanes', value: 'one' } : { op: 'unset', at, key: 'lanes' },
      ];
      const o = spec.patch('Lanes', (_root, revision) => ({ base: revision, ops }), {
        ...(ctx.origin ? { origin: ctx.origin } : {}),
      });
      const done = ctx.lanes === 'one' ? 'One lane' : 'Packed';
      return o instanceof Promise ? o.then(said(done)) : said(done)(o);
    },
  });

  return {
    commands,
    composed(v) {
      viewer = v;
    },
    open: (tracks) => live()?.open(tracks) ?? false,
    dispose() {
      disposed = true;
      stop();
      body?.dispose();
    },
  };
}

/** A section's track groups folded and opened (_spec/README.md "Behaviors", Groups). */
export const GroupsProvider: BehaviorProvider<undefined, GroupsBehavior> = {
  name: 'groups',
  create: ({ spec, view }) => groupsBehavior(spec, view),
};
