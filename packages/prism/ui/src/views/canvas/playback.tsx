/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { Markdown } from '../../components/markdown/markdown.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { For, Show } from '../../dom/control-flow.js';
import { effect } from '../../dom/scope.js';
import { restartable } from '../../dom/scoped.js';
import { keyedSelector } from '../../reactive/keyed-selector.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive, State } from '../../reactive/state.js';
import type { Child } from '../../dom/jsx-types.js';
import type {
  PeriodicBoxMode,
  PeriodicBoxModel,
  SnapshotStepperModel,
  SnapshotTransitionModel,
  TrajectoryControlsModel,
} from '../../model/playback.js';
import { openContextMenu } from '../../components/context-menu/context-menu.js';
import { classes } from './_generated/playback.module.css.js';

/** Up to this many snapshots, the position shows as dots; past it, as "n / N" opening a list of titles. */
const MaxDots = 7;

/**
 * Snapshot playback (the session's story) at the top left of the canvas: the current snapshot's title, its position,
 * previous, next, the slideshow's play (from two snapshots on), then its description. Only when the session has
 * snapshots.
 */
export function SnapshotCard(props: {
  model: SnapshotStepperModel;
  /** In a card already (an overlay's): no glass of its own, so there is one frame. */
  bare?: boolean;
}) {
  return (
    <Show when={derive(props.model.count, (c) => c > 0)}>
      {() => (
        <section class={props.bare ? classes.bare : classes.card} aria-label="Snapshots">
          <Snapshots model={props.model} />
        </section>
      )}
    </Show>
  );
}

/**
 * The canvas's playback controls in one card, a row each: the snapshots (step through them, play), the current
 * snapshot's transition (a scrubber), the trajectory's models. Each row shows only when it has something to play; with
 * none, the card is not there.
 */
export function PlaybackDeck(props: { children?: Child }) {
  return (
    <section class={classes.deck} aria-label="Playback">
      {props.children}
    </section>
  );
}

/** Trajectory playback alone: the deck with its row. */
export function TrajectoryBar(props: { model: TrajectoryControlsModel }) {
  return (
    <PlaybackDeck>
      <TrajectoryRow model={props.model} />
    </PlaybackDeck>
  );
}

/**
 * The trajectory's row: the models of one structure (named, since it is that structure's model index), a scrubber
 * and play; `compact`, only its name, the frame and play (where the full controls are at hand elsewhere: the viewer's
 * Structure panel). Only for a multi-model structure.
 */
export function TrajectoryRow(props: { model: TrajectoryControlsModel; compact?: boolean }) {
  return (
    <Show when={derive(props.model.count, (c) => c > 1)}>
      {() => <Trajectory model={props.model} compact={props.compact} />}
    </Show>
  );
}

/**
 * The snapshots' row (v5's viewport snapshot controls): the position (the current one's title on hover; a list of all
 * of them), previous, play, next, and `actions` (the transition's toggle, mute). From two snapshots on; the Story card
 * has the title and the description.
 */
export function SnapshotRow(props: {
  model: SnapshotStepperModel;
  /** At the row's end, after next (a story's audio: mute). */
  actions?: () => Child;
}) {
  const { model } = props;
  return (
    <Show when={derive(model.count, (c) => c > 1)}>
      {() => {
        const listOpen = observableState(false);
        return (
          <div data-row>
            <div class={classes.bar} role="group" aria-label="Snapshots">
              <span class={classes.lead}>
                <Icon name="snapshots" size={16} />
              </span>
              {/* No title: quiet over the canvas (the Story card has it); the position says it on hover, and opens
                  the list of all of them. */}
              <button
                type="button"
                class={classes.position}
                aria-expanded={listOpen}
                aria-label={derive(model.title, (t) => (t ? `Snapshots: ${t}` : 'Snapshots'))}
                data-tip={derive(model.title, (t) => t ?? 'No snapshot applied')}
                onClick={() => listOpen.set(!listOpen.value)}
              >
                {derive([model.index, model.count], (i, c) => `${i < 0 ? '–' : i + 1} / ${c}`)}
              </button>
              <Button
                variant="ghost"
                size="sm"
                icon={<Icon name="chevron-left" />}
                title="Previous snapshot"
                onClick={() => model.previous()}
              />
              <PlayButton playing={model.isPlaying} title="Play snapshots" onClick={() => model.togglePlay()} />
              <Button
                variant="ghost"
                size="sm"
                icon={<Icon name="chevron-right" />}
                title="Next snapshot"
                onClick={() => model.next()}
              />
              {props.actions?.()}
            </div>
            <Show when={listOpen}>{() => <SnapshotList model={model} onClose={() => listOpen.set(false)} />}</Show>
          </div>
        );
      }}
    </Show>
  );
}

/** A mute button: the speaker, struck through while muted (pressed). */
export function MuteButton(props: { muted: ReadonlyState<boolean>; onClick(): void }) {
  return (
    <Button
      variant="ghost"
      size="sm"
      icon={
        <Show when={props.muted} fallback={() => <Icon name="volume-2" />}>
          {() => <Icon name="volume-x" />}
        </Show>
      }
      title={derive(props.muted, (m) => (m ? 'Unmute' : 'Mute'))}
      pressed={props.muted}
      onClick={() => props.onClick()}
    />
  );
}

/**
 * A story's audio in the playback card when it has no snapshots row to sit in (one snapshot): its name and mute.
 */
export function AudioRow(props: { muted: ReadonlyState<boolean>; onToggleMute(): void }) {
  return (
    <div class={classes.bar} role="group" aria-label="Audio" data-row>
      <span class={classes.lead}>
        <Icon name="volume-2" size={16} />
      </span>
      <span class={classes.rowTitle}>Audio</span>
      <MuteButton muted={props.muted} onClick={props.onToggleMute} />
    </div>
  );
}

/** Seconds with a tenth ("1.2 s"), the transition's time. */
const seconds = (ms: number) => (ms / 1000).toFixed(1);

/**
 * The current snapshot's transition (v5's snapshot animation slider): a scrubber through it, where it is, play from its
 * start or stop. Only when the snapshot has one.
 */
export function TransitionRow(props: { model: SnapshotTransitionModel }) {
  const { model } = props;
  return (
    <Show when={derive(model.duration, (d) => d > 0)}>
      {() => {
        const progress = derive([model.time, model.duration], (t, d) => `${d > 0 ? (100 * Math.min(t, d)) / d : 0}%`);
        return (
          <div class={classes.bar} role="group" aria-label="Transition" data-row>
            <span class={classes.lead}>
              <Icon name="film" size={16} />
            </span>
            <input
              class={classes.scrubber}
              type="range"
              min="0"
              max={derive(model.duration, String)}
              step={derive(model.step, (s) => String(Math.max(1, Math.round(s))))}
              value={derive(model.time, (t) => String(Math.round(t)))}
              style={{ '--msp-progress': progress }}
              aria-label="Transition time"
              onInput={(e) => model.seek(Number(e.currentTarget.value))}
            />
            <span class={classes.time}>
              {derive([model.time, model.duration], (t, d) => `${seconds(t)} / ${seconds(d)} s`)}
            </span>
            <Button
              variant="ghost"
              size="sm"
              icon={
                <Show when={model.isPlaying} fallback={() => <Icon name="play" />}>
                  {() => <Icon name="square" />}
                </Show>
              }
              title={derive(model.isPlaying, (p) => (p ? 'Stop the transition' : 'Play the transition'))}
              onClick={() => model.togglePlay()}
            />
          </div>
        );
      }}
    </Show>
  );
}

function Snapshots(props: { model: SnapshotStepperModel }) {
  const { model } = props;
  const listOpen = observableState(false);
  const many = derive(model.count, (c) => c > MaxDots);
  // The title and the description change together (another snapshot): taken once both have (the task after), so a
  // heading is never compared with the other snapshot's title.
  const description = observableState<string | undefined>(undefined);
  let taken: string | undefined;
  const timer = restartable(() => description.set(taken));
  effect(
    derive([model.description, model.title], (d, t) => (d ? withoutTitle(d, t ?? '') : d)),
    (d) => {
      taken = d;
      timer.start(0);
    },
  );
  return (
    <div class={classes.snapshots}>
      <div class={classes.bar}>
        <span class={classes.title}>{model.title}</span>
        {/* One snapshot (a single-state story) has its title and description only: nothing to step through. */}
        <Show when={derive(model.count, (c) => c > 1)}>
          {() => (
            <>
              <Show when={many} fallback={() => <Dots model={model} />}>
                {() => (
                  <button
                    type="button"
                    class={classes.position}
                    aria-expanded={listOpen}
                    aria-label="Snapshots"
                    onClick={() => listOpen.set(!listOpen.value)}
                  >
                    {derive([model.index, model.count], (i, c) => `${i + 1} / ${c}`)}
                  </button>
                )}
              </Show>
              <Button
                variant="ghost"
                size="sm"
                icon={<Icon name="chevron-left" />}
                title="Previous snapshot"
                onClick={() => model.previous()}
              />
              <Button
                variant="ghost"
                size="sm"
                icon={<Icon name="chevron-right" />}
                title="Next snapshot"
                onClick={() => model.next()}
              />
              <PlayButton playing={model.isPlaying} title="Play snapshots" onClick={() => model.togglePlay()} />
            </>
          )}
        </Show>
      </div>
      <Show when={listOpen}>{() => <SnapshotList model={model} onClose={() => listOpen.set(false)} />}</Show>
      <Show when={description.asReadonly}>
        {(source) => (
          <ScrollArea class={classes.description} contentClass={classes.descriptionContent}>
            <Show
              when={model.plainText ?? false}
              fallback={() => (
                <Markdown source={source} onLink={(href) => model.openLink(href)} host={model.markdown} />
              )}
            >
              {() => <p class={classes.plain}>{source}</p>}
            </Show>
          </ScrollArea>
        )}
      </Show>
    </div>
  );
}

/**
 * A description without its first heading when that says nothing the snapshot's title (shown in the bar) does not: the
 * same words, or fewer of them ("ATP Binding" under "ATP Binding [1/2]"), case, markdown marks, punctuation and a part
 * counter ("[1/2]", "(2/5)") aside. A heading that says more ("The Story of BCR-ABL: A Kinase Out of Control" under "A
 * Kinase Out of Control") stays.
 */
export function withoutTitle(markdown: string, title: string): string {
  const m = /^\s*#{1,6}[ \t]+(.+?)[ \t]*#*[ \t]*(?:\r?\n|$)/.exec(markdown);
  if (!m) return markdown;
  const heading = words(m[1]);
  const said = words(title);
  if (
    !heading ||
    !(
      said === heading ||
      said.startsWith(`${heading} `) ||
      said.includes(` ${heading} `) ||
      said.endsWith(` ${heading}`)
    )
  )
    return markdown;
  return markdown.slice(m[0].length).replace(/^\s+/, '');
}

/** A heading's or a title's words: lower case, no markdown marks, punctuation or part counter. */
function words(text: string): string {
  return text
    .replace(/[[(]\s*\d+\s*\/\s*\d+\s*[\])]/g, ' ')
    .replace(/[*_`~]/g, '')
    .toLowerCase()
    .replace(/[^\p{L}\p{N}]+/gu, ' ')
    .trim();
}

/** One dot per snapshot, the current one filled; each jumps to its snapshot (its title as the tooltip). */
function Dots(props: { model: SnapshotStepperModel }) {
  const { model } = props;
  return (
    <span class={classes.dots} role="tablist" aria-label="Snapshots">
      <For each={model.titles}>
        {(title, index) => (
          <button
            type="button"
            role="tab"
            class={classes.dot}
            data-tip={title}
            aria-label={title}
            aria-selected={derive([model.index, index], (i, j) => i === j)}
            onClick={() => model.goTo(index.value)}
          />
        )}
      </For>
    </span>
  );
}

/** Many snapshots: their titles as a list under the head, filterable when long; the current one marked. */
function SnapshotList(props: { model: SnapshotStepperModel; onClose(): void }) {
  const { model } = props;
  const filter = observableState('');
  const rows = derive([model.titles, filter.asReadonly], (titles, f) => {
    const q = f.trim().toLowerCase();
    return titles.map((title, index) => ({ title, index })).filter((r) => !q || r.title.toLowerCase().includes(q));
  });
  const current = keyedSelector(model.index);
  return (
    <div class={classes.list}>
      <Show when={derive(model.count, (c) => c > 12)}>
        {() => (
          <TextInput
            search
            size="sm"
            aria-label="Filter snapshots"
            placeholder="Filter"
            value={filter}
            onInput={(v) => filter.set(v)}
          />
        )}
      </Show>
      <ScrollArea class={classes.listRows} contentClass={classes.listRowsContent} role="listbox" aria-label="Snapshots">
        {/* The current one, told to the two rows it moves between (view-models.md §3). */}
        <For each={rows} by={(r) => r.index}>
          {(r) => (
            <button
              type="button"
              role="option"
              class={classes.listRow}
              aria-selected={current.key(r.index)}
              onClick={() => {
                model.goTo(r.index);
                props.onClose();
              }}
            >
              <span class={classes.listIndex}>{String(r.index + 1)}</span>
              <span>{r.title}</span>
            </button>
          )}
        </For>
      </ScrollArea>
    </div>
  );
}

const ModelsUnit: ReadonlyState<'models' | 'frames'> = new State<'models' | 'frames'>('models').asReadonly;

function Trajectory(props: { model: TrajectoryControlsModel; compact?: boolean }) {
  const { model } = props;
  // An ensemble's models or a simulation's frames: what the controls say (it follows the structure the bar follows).
  const unitOf = model.unit ?? ModelsUnit;
  const Unit = derive(unitOf, (u) => (u === 'models' ? 'Models' : 'Frames'));
  const one = derive(unitOf, (u) => (u === 'models' ? 'model' : 'frame'));
  const shownLabel = derive(Unit, (u) => `${u.slice(0, -1)} shown`);
  if (props.compact) {
    return (
      <div class={classes.bar} role="group" aria-label={Unit} data-row>
        <span class={classes.lead}>
          <Icon name="trajectory" size={16} />
        </span>
        <span class={classes.trajectoryLabel}>{model.label}</span>
        <span class={classes.counter}>{derive([model.index, model.count], (i, c) => `${i + 1} / ${c}`)}</span>
        <PlayButton
          playing={model.isPlaying}
          title={derive(unitOf, (u) => `Play ${u}`)}
          onClick={() => model.togglePlay()}
        />
      </div>
    );
  }
  // One model back or forth, past the end to the start (as v5's arrows).
  const step = (by: number) => {
    const c = model.count.value;
    if (c > 1) model.seek((((model.index.value + by) % c) + c) % c);
  };
  // The played part of the track, for the scrubber's fill.
  const progress = derive([model.index, model.count], (i, c) => `${c > 1 ? (100 * i) / (c - 1) : 0}%`);
  return (
    <div class={classes.bar} role="group" aria-label={Unit} data-row>
      <span class={classes.lead}>
        <Icon name="trajectory" size={16} />
      </span>
      <span class={classes.trajectoryLabel}>{model.label}</span>
      <span class={classes.counter}>{derive([model.index, model.count], (i, c) => `${i + 1} / ${c}`)}</span>
      <input
        class={classes.scrubber}
        type="range"
        min="0"
        max={derive(model.count, (c) => String(c - 1))}
        value={derive(model.index, String)}
        style={{ '--msp-progress': progress }}
        aria-label={shownLabel}
        // Said as people count them ("Model 3 of 10"), not the 0-based value.
        aria-valuetext={derive([Unit, model.index, model.count], (u, i, c) => `${u.slice(0, -1)} ${i + 1} of ${c}`)}
        onInput={(e) => model.seek(Number(e.currentTarget.value))}
      />
      <Button
        variant="ghost"
        size="sm"
        icon={<Icon name="chevron-left" />}
        title={derive(one, (o) => `Previous ${o}`)}
        onClick={() => step(-1)}
      />
      <PlayButton
        playing={model.isPlaying}
        title={derive(unitOf, (u) => `Play ${u}`)}
        onClick={() => model.togglePlay()}
      />
      <Button
        variant="ghost"
        size="sm"
        icon={<Icon name="chevron-right" />}
        title={derive(one, (o) => `Next ${o}`)}
        onClick={() => step(1)}
      />
      {model.periodic ? <PeriodicBoxButton model={model.periodic} /> : null}
    </div>
  );
}

const PeriodicBoxLabels: Record<PeriodicBoxMode, string> = {
  whole: 'Molecules whole',
  box: 'Molecules in the box',
  'box-centered': 'In the box, polymers at the center',
  nojump: 'No jumps between frames',
  off: 'As in the file',
};

/** The periodic box's menu: shown only for a structure in a simulation box; the mode now is marked. */
function PeriodicBoxButton(props: { model: PeriodicBoxModel }) {
  const { model } = props;
  let close: (() => void) | undefined;
  return (
    <Show when={derive(model.mode, (m) => m !== undefined)}>
      {() => {
        const button = (
          <Button
            variant="ghost"
            size="sm"
            icon={<Icon name="periodic-box" />}
            title={derive(model.mode, (m) => `Periodic box: ${m ? PeriodicBoxLabels[m] : ''}`)}
            onClick={(e: MouseEvent) => {
              const r = (e.currentTarget as HTMLElement).getBoundingClientRect();
              close?.();
              const now = model.mode.value;
              close = openContextMenu(
                (e.currentTarget as HTMLElement).closest<HTMLElement>('.msp-root') ?? document.body,
                { x: r.left, y: r.bottom },
                [
                  ...(Object.keys(PeriodicBoxLabels) as PeriodicBoxMode[]).map((m) => ({
                    id: m,
                    label: PeriodicBoxLabels[m],
                    icon: m === now ? ('check' as const) : undefined,
                    run: () => model.set(m),
                  })),
                  // the neighbouring boxes: a switch, after the modes
                  ...(model.images
                    ? [
                        {
                          id: 'images',
                          label: 'Periodic images',
                          icon: model.images.shown.value ? ('check' as const) : undefined,
                          run: () => model.images!.set(!model.images!.shown.value),
                        },
                      ]
                    : []),
                ],
              );
            }}
          />
        );
        return button;
      }}
    </Show>
  );
}

function PlayButton(props: {
  playing: ReadonlyState<boolean>;
  title: string | ReadonlyState<string>;
  onClick(): void;
}) {
  const title = typeof props.title === 'string' ? new State(props.title).asReadonly : props.title;
  return (
    <Button
      variant="ghost"
      size="sm"
      icon={
        <Show when={props.playing} fallback={() => <Icon name="play" />}>
          {() => <Icon name="pause" />}
        </Show>
      }
      title={derive([props.playing, title], (p, t) => (p ? 'Pause' : t))}
      pressed={props.playing}
      onClick={() => props.onClick()}
    />
  );
}
