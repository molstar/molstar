/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { holdCanvas } from '../../sources/canvas-presence.js';
import { OptionalParam } from '../../api/params.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { PluginState } from '@molstar/plugin/state';
import { AnimateStateSnapshotTransition } from '@molstar/plugin/state/animation/built-in/state-snapshots';
import { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import { type ApiContext, ApiArea, ApiResult, type ApiFailure } from '../../api/area.js';
import {
  createMarkdownSource,
  createSnapshotsSource,
  type MarkdownSource,
  type SnapshotsSource,
  type SnapshotsState,
  type SnapshotTransitionState,
} from './source.js';
import type { ReadonlyState } from '@molstar/prism-ui/core';

// The `snapshot` area (api.md §2; viewer-spec.md §3.1, the `snapshots` feature): the plugin's snapshot manager
// as operations. Its checkpoints never use the manager (sequence.md §7.2) and the list is not part of Undo ("what it
// never touches: … the Snapshots list"), so changes to the list are `undoable: false` (no checkpoint, recorded as
// runs). Showing one (apply, next, previous) is no step either (David, 2026-10-20): `undoable: false`, so it runs in
// order with the steps (Undo is refused while one is applied, a step never runs against a half-applied scene), takes no checkpoint, and
// leaves the steps before it undoable; Undo brings back the snapshot that was current with the scene it restores
// (sequence.md §7.3).

type Entry = PluginStateSnapshotManager.Entry;

const Person = { kind: 'person' } as const;

/** One entry as `snapshot.list` reports it. */
export interface SnapshotListEntry {
  readonly id: string;
  /** 1-based, as `snapshot.move` takes it. */
  readonly position: number;
  readonly name: string;
  readonly key?: string;
  readonly description?: string;
  readonly current: boolean;
  /** How long playback shows it (ms). */
  readonly showFor?: number;
}

export interface SnapshotList {
  readonly entries: readonly SnapshotListEntry[];
  readonly playing: boolean;
}

const EntryParam = (description: string) => PD.Text('', { description });
const ByIdNameOrKey = 'A snapshot by id, name or key';

const Options = () =>
  PD.Value<Record<string, unknown>>(
    {},
    {
      isHidden: true,
      description: "What the snapshot includes, as Mol*'s snapshot options: data, camera, canvas3d, interactivity, …",
    },
  );

/** How long playback shows a snapshot; absent: `absent` (api.md §1). */
const ShowFor = (absent = 'the default') =>
  OptionalParam(
    PD.Numeric(
      1500,
      { min: 100, max: 60000, step: 100 },
      { description: `How long playback shows it, in ms (absent: ${absent})` },
    ),
  );

const AddParams = {
  name: PD.Text('', { description: 'Its name; empty: Snapshot <n>' }),
  description: PD.Text('', { description: 'Markdown shown with it; link another with [title](#key)' }),
  key: PD.Text('', { description: 'A key to link it from descriptions; unique' }),
  showFor: ShowFor(),
  options: Options(),
};

const UpdateParams = {
  entry: EntryParam(`${ByIdNameOrKey}; empty: the current one`),
  // Each given is set (empty text clears it), each absent kept (api.md §1).
  name: OptionalParam(PD.Text('', { description: 'A new name ("": none)' })),
  description: OptionalParam(PD.Text('', { description: 'A new description in markdown ("": none)' })),
  key: OptionalParam(PD.Text('', { description: 'A new key ("": none)' })),
  showFor: OptionalParam(
    PD.Numeric(1500, { min: 100, max: 60000, step: 100 }, { description: 'How long playback shows it, in ms' }),
  ),
};

const ReplaceParams = {
  entry: EntryParam(`${ByIdNameOrKey}; empty: the current one`),
  showFor: ShowFor('keep its own'),
  options: Options(),
};

const MoveParams = {
  entry: EntryParam(ByIdNameOrKey),
  position: PD.Numeric(1, { min: 1, step: 1 }, { description: 'Where it goes in the list (1 is the first)' }),
};

const PlayParams = {
  delay: OptionalParam(
    PD.Numeric(
      1500,
      { min: 100, max: 60000, step: 100 },
      { description: "How long to show a snapshot that has no time of its own, in ms (absent: the manager's)" },
    ),
  ),
};

const SeekParams = {
  time: PD.Numeric(
    0,
    { min: 0, max: 600000, step: 1 },
    { description: "Where in the current snapshot's transition, in ms" },
  ),
};

const TransitionParams = {
  play: PD.Boolean(true, { description: "Play the current snapshot's transition from its start; false: stop it" }),
};

const LinkParams = {
  href: PD.Text('', {
    description: "The link: a Mol* markdown command, '!name=value&…' (focus-refs, highlight-refs, …)",
  }),
  event: PD.Select<'click' | 'enter' | 'leave'>('click', [
    ['click', 'Click'],
    ['enter', 'Pointer enters'],
    ['leave', 'Pointer leaves'],
  ]),
};

export type AddSnapshotParams = typeof AddParams;
export type UpdateSnapshotParams = typeof UpdateParams;
export type ReplaceSnapshotParams = typeof ReplaceParams;
export type MoveSnapshotParams = typeof MoveParams;
export type PlayParams = typeof PlayParams;

const labelOf = (e: Entry | undefined, i: number) => `'${e?.name || `Snapshot ${i + 1}`}'`;
const entriesOf = (ctx: ApiContext) => ctx.plugin.managers.snapshot.state.entries.toArray();
const none = () => ApiResult.fail('not-found', 'There are no snapshots.', 'Add one first (snapshot.add).');

/** An entry by id, then by name (any case), then by key; empty means the current one when `orCurrent`. */
function find(ctx: ApiContext, ref: string, orCurrent: boolean): { entry: Entry; index: number } | ApiFailure {
  const entries = entriesOf(ctx);
  if (!entries.length) return none();
  const text = ref.trim();
  if (!text) {
    const id = ctx.plugin.managers.snapshot.state.current;
    const index = entries.findIndex((e) => e.snapshot.id === id);
    if (orCurrent && index >= 0) return { entry: entries[index], index };
    return ApiResult.fail(
      'invalid',
      orCurrent ? 'No snapshot is current.' : 'Which snapshot?',
      'Give its id, name or key (snapshot.list).',
    );
  }
  let index = entries.findIndex((e) => e.snapshot.id === text);
  if (index < 0) {
    const lower = text.toLowerCase();
    const named = entries.flatMap((e, i) => (e.name?.toLowerCase() === lower ? [i] : []));
    if (named.length > 1) {
      return ApiResult.fail(
        'ambiguous',
        `${named.length} snapshots are named '${text}'.`,
        `Use an id: ${named.map((i) => `${i + 1}. ${entries[i].snapshot.id}`).join('; ')}.`,
      );
    }
    index = named[0] ?? entries.findIndex((e) => e.key === text);
  }
  if (index < 0) return ApiResult.fail('not-found', `There is no snapshot '${text}'.`, 'List them with snapshot.list.');
  return { entry: entries[index], index };
}

const isFailure = (x: unknown): x is ApiFailure => ApiResult.isFailure(x);

/** A key is unique among the snapshots. */
function keyTaken(ctx: ApiContext, key: string, except?: Entry): boolean {
  return !!key && entriesOf(ctx).some((e) => e !== except && e.key === key);
}

/** The state of a new or replacing snapshot: what Mol*'s Add and Replace take, without the canvas where there is none. */
async function capture(ctx: ApiContext, options: Record<string, unknown>, showFor: number | undefined) {
  const plugin = ctx.plugin;
  const params: PluginState.SnapshotParams = { ...options };
  if (showFor !== undefined) params.durationInMs = showFor;
  // `getSnapshot` reads the camera from the canvas; a thumbnail needs one.
  if (!plugin.canvas3d) {
    params.camera = false;
    params.image = false;
  }
  const snapshot = plugin.state.getSnapshot(params);
  const wantsImage = params.image ?? plugin.state.snapshotParams.value.image;
  let image: Awaited<ReturnType<typeof PluginStateSnapshotManager.getCanvasImageAsset>> | undefined;
  if (wantsImage) {
    // The thumbnail is rendered on purpose: the canvas is held while it is (README decision 29).
    // A context the browser will not give back leaves the snapshot without its thumbnail.
    const unhold = await holdCanvas(plugin, { signal: ctx.signal }).catch(() => undefined);
    if (unhold) {
      try {
        image = await PluginStateSnapshotManager.getCanvasImageAsset(plugin, `${snapshot.id}-image.png`);
      } finally {
        unhold();
      }
    }
  }
  return { snapshot, image };
}

/** The manager keeps `isPlaying` and its default delay in protected state (its own TODO: not editable). */
function setManagerState(
  ctx: ApiContext,
  patch: { current?: string; entries?: PluginStateSnapshotManager['state']['entries']; nextSnapshotDelayInMs?: number },
) {
  (ctx.plugin.managers.snapshot as unknown as { updateState(s: object): boolean }).updateState(patch);
}

export class SnapshotApi extends ApiArea {
  readonly id = 'snapshot';

  #source: SnapshotsSource | undefined;
  #markdown: MarkdownSource | undefined;

  /**
   * The snapshots as a state the views read (their list, the current one, whether it plays); made on first use and
   * released with the area. Not an operation: a read of the plugin that stays in this feature's source.
   */
  get state(): ReadonlyState<SnapshotsState> {
    return (this.#source ??= createSnapshotsSource(this.ctx.plugin)).state;
  }

  /** The current snapshot's transition as the canvas scrubber reads it (its length, where it is, whether it plays). */
  get transitionState(): ReadonlyState<SnapshotTransitionState> {
    return (this.#source ??= createSnapshotsSource(this.ctx.plugin)).transition;
  }

  /**
   * The plugin's markdown extensions as the descriptions read them (which links are commands, images, the audio
   * player); made on first use. Running a command is `snapshot.link`.
   */
  get markdown(): MarkdownSource {
    return (this.#markdown ??= createMarkdownSource(this.ctx.plugin));
  }

  /**
   * A story's own navigation (a description's `!apply-snapshot` / `!next-snapshot`, a clickable primitive's
   * `snapshot_key`, the keyboard's next and previous) goes through these operations as the person, in order with the
   * history like the panel's (sequence.md §7.3): the plugin's manager hands it here.
   */
  readonly #unnavigate = this.ctx.plugin.managers.snapshot.setNavigator({
    applyKey: (key) => {
      const id = this.ctx.plugin.managers.snapshot.state.entries.find((e) => e.key === key)?.snapshot.id;
      return id ? this.ctx.operations.run(this.apply.definition, { entry: id }, { origin: Person }) : undefined;
    },
    applyNext: (dir) =>
      this.ctx.operations.run((dir > 0 ? this.next : this.previous).definition, {}, { origin: Person }),
  });

  dispose(): void {
    this.#unnavigate();
    this.#source?.dispose();
    this.#source = undefined;
    this.#markdown?.dispose();
    this.#markdown = undefined;
    // A running slideshow holds a timer of its own.
    try {
      const manager = this.ctx.plugin.managers.snapshot;
      if (manager.state.isPlaying) manager.stop().catch(() => {});
    } catch {
      // The plugin is gone already.
    }
  }

  readonly list = this.operation({
    id: 'snapshot.list',
    title: 'List the snapshots',
    params: {},
    read: true,
    async run(_, ctx) {
      const { entries, current, isPlaying } = ctx.plugin.managers.snapshot.state;
      const list: SnapshotList = {
        playing: isPlaying,
        entries: entries.toArray().map((e, i) => ({
          id: e.snapshot.id,
          position: i + 1,
          name: e.name ?? '',
          ...(e.key ? { key: e.key } : {}),
          ...(e.description ? { description: e.description } : {}),
          current: e.snapshot.id === current,
          ...(e.snapshot.durationInMs !== undefined ? { showFor: e.snapshot.durationInMs } : {}),
        })),
      };
      if (!list.entries.length) return ApiResult.ok('There are no snapshots.', {}, { value: list });
      return ApiResult.ok(
        list.entries
          .map((e) => `${e.position}. ${e.name || `Snapshot ${e.position}`}${e.current ? ' (current)' : ''} [${e.id}]`)
          .join('\n'),
        {},
        { value: list },
      );
    },
  });

  readonly add = this.operation({
    id: 'snapshot.add',
    title: 'Add a snapshot',
    params: AddParams,
    undoable: false,
    async run({ name, description, key, showFor, options }, ctx) {
      if (keyTaken(ctx, key)) return ApiResult.fail('invalid', `The key '${key}' is taken.`, 'A key is unique.');
      const manager = ctx.plugin.managers.snapshot;
      const title = name.trim() || `Snapshot ${manager.state.entries.size + 1}`;
      const { snapshot, image } = await capture(ctx, options, showFor);
      manager.add(
        PluginStateSnapshotManager.Entry(snapshot, {
          name: title,
          ...(description.trim() ? { description: description.trim() } : {}),
          ...(key.trim() ? { key: key.trim() } : {}),
          ...(image ? { image } : {}),
        }),
      );
      return ApiResult.ok(
        `Added snapshot '${title}' (${manager.state.entries.size} in all).`,
        {},
        { value: { id: snapshot.id, name: title, position: manager.state.entries.size } },
      );
    },
  });

  readonly apply = this.operation({
    id: 'snapshot.apply',
    // Stepping through a story: the card and the playback row show it, no message (David, 2026-10-15).
    quiet: true,
    title: 'Apply a snapshot',
    params: { entry: EntryParam(ByIdNameOrKey) },
    undoable: false,
    async run({ entry }, ctx) {
      const found = find(ctx, entry, false);
      if (isFailure(found)) return found;
      const manager = ctx.plugin.managers.snapshot;
      const snapshot = manager.getEntry(found.entry.snapshot.id)?.snapshot;
      if (!snapshot) return ApiResult.fail('not-found', `There is no snapshot '${entry}'.`);
      await ctx.plugin.state.setSnapshot(snapshot);
      // Current only once it is applied: a failed apply leaves the list's current entry as it was.
      manager.setCurrent(snapshot.id);
      return ApiResult.ok(
        `Applied snapshot ${labelOf(found.entry, found.index)}.`,
        { view: true, sequence: true },
        { value: { id: snapshot.id, position: found.index + 1 } },
      );
    },
  });

  readonly update = this.operation({
    id: 'snapshot.update',
    title: 'Edit a snapshot',
    surface: { kind: 'agent-only', reason: "The Snapshots tab's rows." },
    params: UpdateParams,
    undoable: false,
    async run(v, ctx) {
      const { entry } = v;
      const { name, description, key, showFor } = v as {
        name?: string;
        description?: string;
        key?: string;
        showFor?: number;
      };
      const found = find(ctx, entry, true);
      if (isFailure(found)) return found;
      const e = found.entry;
      const text = (given: string | undefined, current: string | undefined) =>
        given === undefined ? current : given.trim() || undefined;
      const next = {
        name: text(name, e.name),
        description: text(description, e.description),
        key: text(key, e.key),
      };
      const time = showFor ?? e.snapshot.durationInMs;
      if (
        next.name === e.name &&
        next.description === e.description &&
        next.key === e.key &&
        time === e.snapshot.durationInMs
      ) {
        return ApiResult.already(`Snapshot ${labelOf(e, found.index)} is as it was.`);
      }
      if (next.key && keyTaken(ctx, next.key, e)) {
        return ApiResult.fail('invalid', `The key '${next.key}' is taken.`, 'A key is unique.');
      }
      // The snapshot's own playback time is part of its state, which the entry holds by reference.
      if (time !== undefined) e.snapshot.durationInMs = time;
      ctx.plugin.managers.snapshot.update(e, { ...next, descriptionFormat: e.descriptionFormat });
      return ApiResult.ok(`Updated snapshot ${labelOf({ ...e, name: next.name }, found.index)}.`, {});
    },
  });

  readonly replace = this.operation({
    id: 'snapshot.replace',
    title: 'Replace a snapshot with the current view',
    surface: { kind: 'agent-only', reason: "The Snapshots tab's rows." },
    params: ReplaceParams,
    undoable: false,
    async run({ entry, showFor, options }, ctx) {
      const found = find(ctx, entry, true);
      if (isFailure(found)) return found;
      // The time it is shown for stays the entry's own (the options' time is for new snapshots).
      const time = (showFor as number | undefined) ?? found.entry.snapshot.durationInMs ?? 0;
      const { snapshot, image } = await capture(ctx, options, time);
      // The id changes with the state it holds (the manager's rule).
      ctx.plugin.managers.snapshot.replace(found.entry.snapshot.id, snapshot, { image });
      return ApiResult.ok(
        `Replaced snapshot ${labelOf(found.entry, found.index)} with the current state.`,
        {},
        { value: { id: snapshot.id, position: found.index + 1 } },
      );
    },
  });

  readonly move = this.operation({
    id: 'snapshot.move',
    title: 'Move a snapshot',
    surface: { kind: 'agent-only', reason: "The Snapshots tab's rows." },
    params: MoveParams,
    undoable: false,
    async run({ entry, position }, ctx) {
      const found = find(ctx, entry, false);
      if (isFailure(found)) return found;
      const count = entriesOf(ctx).length;
      const to = Math.round(position) - 1;
      if (to < 0 || to >= count) {
        return ApiResult.fail('invalid', `Position ${position} is not in the list.`, `There are ${count}.`);
      }
      if (to === found.index) {
        return ApiResult.already(`Snapshot ${labelOf(found.entry, found.index)} is already at ${position}.`);
      }
      const manager = ctx.plugin.managers.snapshot;
      // In one step: the manager's `move` swaps neighbors (a search of the list and a change told, each), so a move
      // across the list was its length × that. The current one stays (the one applied last: moving does not apply).
      const entries = manager.state.entries.delete(found.index).insert(to, found.entry);
      setManagerState(ctx, { current: manager.state.current, entries });
      manager.events.changed.next(void 0);
      return ApiResult.ok(`Moved snapshot ${labelOf(found.entry, found.index)} to ${position}.`, {});
    },
  });

  readonly remove = this.operation({
    id: 'snapshot.remove',
    title: 'Remove a snapshot',
    surface: { kind: 'agent-only', reason: "The Snapshots tab's rows." },
    params: { entry: EntryParam(ByIdNameOrKey) },
    undoable: false,
    approval: 'destructive',
    async run({ entry }, ctx) {
      const found = find(ctx, entry, false);
      if (isFailure(found)) return found;
      ctx.plugin.managers.snapshot.remove(found.entry.snapshot.id);
      return ApiResult.ok(
        `Removed snapshot ${labelOf(found.entry, found.index)} (${entriesOf(ctx).length} left).`,
        {},
        { value: { id: found.entry.snapshot.id } },
      );
    },
  });

  readonly clear = this.operation({
    id: 'snapshot.clear',
    title: 'Remove all snapshots',
    surface: { kind: 'agent-only', reason: "The Snapshots tab's rows." },
    params: {},
    undoable: false,
    approval: 'destructive',
    async run(_, ctx) {
      const count = entriesOf(ctx).length;
      if (!count) return ApiResult.already('There are no snapshots.');
      ctx.plugin.managers.snapshot.clear();
      return ApiResult.ok(`Removed ${count} snapshot${count === 1 ? '' : 's'}.`, {});
    },
  });

  readonly next = this.operation({
    id: 'snapshot.next',
    // Stepping through a story: the card and the playback row show it, no message (David, 2026-10-15).
    quiet: true,
    title: 'Next snapshot',
    params: {},
    undoable: false,
    // Two or more to step through (the arrows' commands are quiet where they cannot run).
    when: (ctx) => (ctx.plugin.managers.snapshot.state.entries.size > 1 ? true : 'There is no story to step through.'),
    run: (_, ctx) => step(ctx, 1),
  });

  readonly previous = this.operation({
    id: 'snapshot.previous',
    // Stepping through a story: the card and the playback row show it, no message (David, 2026-10-15).
    quiet: true,
    title: 'Previous snapshot',
    params: {},
    undoable: false,
    // Two or more to step through (the arrows' commands are quiet where they cannot run).
    when: (ctx) => (ctx.plugin.managers.snapshot.state.entries.size > 1 ? true : 'There is no story to step through.'),
    run: (_, ctx) => step(ctx, -1),
  });

  readonly play = this.operation({
    id: 'snapshot.play',
    title: 'Play the snapshots',
    params: PlayParams,
    history: false,
    async run({ delay }, ctx) {
      const manager = ctx.plugin.managers.snapshot;
      const count = manager.state.entries.size;
      if (!count) return none();
      if (count < 2) return ApiResult.fail('blocked', 'Playing needs at least two snapshots.', 'Add another.');
      if (manager.state.isPlaying) return ApiResult.already('The snapshots are playing.');
      if (delay !== undefined) {
        // The manager keeps the delay for good (and saves it in a session): this play's only.
        const was = manager.state.nextSnapshotDelayInMs;
        setManagerState(ctx, { nextSnapshotDelayInMs: delay });
        const stopped = manager.events.changed.subscribe(() => {
          if (manager.state.isPlaying) return;
          stopped.unsubscribe();
          setManagerState(ctx, { nextSnapshotDelayInMs: was });
        });
      }
      await manager.play();
      // The manager applies each snapshot on its own timer: a presentation, not steps; the history leaves it alone
      // (sources/checkpoint.ts) and Undo after it restores as usual.
      return ApiResult.ok(`Playing ${count} snapshots.`, { view: true });
    },
  });

  /**
   * A command link in a description (v5's markdown commands: camera, focus and highlight refs, query, snapshots, audio),
   * run by the plugin's markdown extensions as v5 runs them. A story's links present it, as playback does: view
   * state, no history step.
   */
  readonly link = this.operation({
    id: 'snapshot.link',
    title: 'Follow a link in a description',
    surface: { kind: 'ui-only', reason: "A description's command links: the person clicks or points at them." },
    params: LinkParams,
    history: false,
    // Pointer-rate (enter and leave): the commands a link runs are recorded as themselves.
    record: false,
    async run({ href, event }, ctx) {
      const extensions = ctx.plugin.managers.markdownExtensions;
      const args = extensions.parseArgs(href);
      if (!args) {
        return ApiResult.fail('invalid', `'${href}' is not a command link.`, "A command link starts with '!'.");
      }
      extensions.tryExecute(event === 'click' ? 'click' : event === 'enter' ? 'mouse-enter' : 'mouse-leave', args);
      const names = Object.keys(args).join(', ');
      return ApiResult.ok(
        event === 'click' ? `Ran ${names}.` : `Ran ${names} (the pointer ${event === 'enter' ? 'entered' : 'left'}).`,
        {
          view: true,
        },
      );
    },
  });

  readonly stop = this.operation({
    id: 'snapshot.stop',
    title: 'Stop playing the snapshots',
    params: {},
    history: false,
    async run(_, ctx) {
      const manager = ctx.plugin.managers.snapshot;
      if (!manager.state.isPlaying) return ApiResult.already('The snapshots are not playing.');
      await manager.stop();
      return ApiResult.ok('Stopped playing the snapshots.', {});
    },
  });

  /**
   * A frame of the current snapshot's transition (Mol*'s `transition`: frames after it is applied), as v5's slider
   * shows it. View state: no step, as playing it is not.
   */
  readonly seek = this.operation({
    id: 'snapshot.seek',
    title: "Show a frame of the snapshot's transition",
    params: SeekParams,
    surface: { kind: 'agent-only', reason: "The playback card's transition scrubber." },
    history: false,
    // A scrubber's steps: not a record entry each.
    record: false,
    when: (ctx) =>
      ctx.plugin.managers.snapshot.current?.snapshot.transition ? true : 'The snapshot has no transition.',
    async run({ time }, ctx) {
      const snapshot = ctx.plugin.managers.snapshot.current?.snapshot;
      if (!snapshot?.transition)
        return ApiResult.fail('blocked', 'The snapshot has no transition.', 'Apply one that has.');
      const animation = ctx.plugin.managers.animation;
      if (animation.isAnimatingStateTransition) await animation.stop();
      const duration = PluginState.getStateTransitionDuration(snapshot) ?? 0;
      const at = Math.max(0, Math.min(time, duration));
      ctx.plugin.managers.snapshot.setSnapshotAnimationFrame(at, true);
      return ApiResult.ok(`At ${Math.round(at)} of ${Math.round(duration)} ms of the transition.`, { view: true });
    },
  });

  /** Plays the current snapshot's transition from its start (v5's replay), or stops it. */
  readonly transition = this.operation({
    id: 'snapshot.transition',
    title: "Play the snapshot's transition",
    params: TransitionParams,
    surface: { kind: 'agent-only', reason: "The playback card's transition row." },
    history: false,
    // The plugin ticks its animations on the canvas's frames.
    modules: ['canvas'],
    when: (ctx) =>
      ctx.plugin.managers.snapshot.current?.snapshot.transition ? true : 'The snapshot has no transition.',
    async run({ play }, ctx) {
      const animation = ctx.plugin.managers.animation;
      if (!play) {
        if (!animation.isAnimatingStateTransition) return ApiResult.already('The transition is not playing.');
        await animation.stop();
        return ApiResult.ok('Stopped the transition.', {});
      }
      if (!ctx.plugin.managers.snapshot.current?.snapshot.transition) {
        return ApiResult.fail('blocked', 'The snapshot has no transition.', 'Apply one that has.');
      }
      if (animation.isAnimatingStateTransition) await animation.stop();
      await animation.play(AnimateStateSnapshotTransition, {});
      return ApiResult.ok("Playing the snapshot's transition.", { view: true });
    },
  });
}

/** The neighbor of the current one, wrapping around as the manager does (none current: the first or the last). */
async function step(ctx: ApiContext, dir: -1 | 1): Promise<ApiResult> {
  const manager = ctx.plugin.managers.snapshot;
  if (!manager.state.entries.size) return none();
  const id = manager.getNextId(manager.state.current, dir);
  const snapshot = id ? manager.getEntry(id)?.snapshot : undefined;
  if (!id || !snapshot) return ApiResult.fail('not-found', 'There is no such snapshot.');
  const index = entriesOf(ctx).findIndex((e) => e.snapshot.id === id);
  await ctx.plugin.state.setSnapshot(snapshot);
  // Current only once it is applied (see `snapshot.apply`).
  manager.setCurrent(id);
  return ApiResult.ok(
    `Applied snapshot ${labelOf(manager.getEntry(id), index)} (${index + 1} of ${manager.state.entries.size}).`,
    { view: true, sequence: true },
    { value: { id, position: index + 1 } },
  );
}
