/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { Markdown } from '../../components/markdown/markdown.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { For, Show } from '../../dom/control-flow.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import type { SnapshotStepperModel, TrajectoryControlsModel } from '../../model/playback.js';
import { classes } from './_generated/playback.module.css.js';

/** Up to this many snapshots, the position shows as dots; past it, as "n / N" opening a list of titles. */
const MaxDots = 7;

/**
 * Snapshot playback (the session's story) at the top left of the canvas: the current snapshot's title, its position,
 * previous, next, the slideshow's play (from two snapshots on), then its description. Only when the session has
 * snapshots.
 */
export function SnapshotCard(props: { model: SnapshotStepperModel }) {
  return (
    <Show when={props.model.count.pipe(map((c) => c > 0))}>
      {() => (
        <section class={classes.card} aria-label="Snapshots">
          <Snapshots model={props.model} />
        </section>
      )}
    </Show>
  );
}

/**
 * Trajectory playback, a separate bar: the models of one structure (named, since it is that structure's model index),
 * a scrubber and play. Only for a multi-model structure.
 */
export function TrajectoryBar(props: { model: TrajectoryControlsModel }) {
  return (
    <Show when={props.model.count.pipe(map((c) => c > 1))}>
      {() => (
        <section class={classes.card} aria-label="Trajectory">
          <Trajectory model={props.model} />
        </section>
      )}
    </Show>
  );
}

function Snapshots(props: { model: SnapshotStepperModel }) {
  const { model } = props;
  const listOpen = observableState(false);
  const many = model.count.pipe(
    map((c) => c > MaxDots),
    distinctUntilChanged(),
  );
  return (
    <div>
      <div class={classes.bar}>
        <span class={classes.title}>{model.title}</span>
        {/* One snapshot (a single-state story) has its title and description only: nothing to step through. */}
        <Show when={model.count.pipe(map((c) => c > 1))}>
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
                    {combineLatest([model.index, model.count]).pipe(map(([i, c]) => `${i + 1} / ${c}`))}
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
      <Show when={model.description}>
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
            aria-selected={combineLatest([model.index, index]).pipe(map(([i, j]) => i === j))}
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
  const rows = combineLatest([model.titles, filter.asReadonly]).pipe(
    map(([titles, f]) => {
      const q = f.trim().toLowerCase();
      return titles.map((title, index) => ({ title, index })).filter((r) => !q || r.title.toLowerCase().includes(q));
    }),
  );
  return (
    <div class={classes.list}>
      <Show when={model.count.pipe(map((c) => c > 12))}>
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
        <For each={rows} key={(r) => r.index}>
          {(r) => (
            <button
              type="button"
              role="option"
              class={classes.listRow}
              aria-selected={model.index.pipe(map((i) => i === r.index))}
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

function Trajectory(props: { model: TrajectoryControlsModel }) {
  const { model } = props;
  // One model back or forth, past the end to the start (as v5's arrows).
  const step = (by: number) => {
    const c = model.count.value;
    if (c > 1) model.seek((((model.index.value + by) % c) + c) % c);
  };
  // The played part of the track, for the scrubber's fill.
  const progress = combineLatest([model.index, model.count]).pipe(
    map(([i, c]) => `${c > 1 ? (100 * i) / (c - 1) : 0}%`),
  );
  return (
    <div class={classes.bar} role="group" aria-label="Models">
      <span class={classes.lead}>
        <Icon name="trajectory" size={16} />
      </span>
      <span class={classes.trajectoryLabel}>{model.label}</span>
      <span class={classes.counter}>
        {combineLatest([model.index, model.count]).pipe(map(([i, c]) => `${i + 1} / ${c}`))}
      </span>
      <input
        class={classes.scrubber}
        type="range"
        min="0"
        max={model.count.pipe(map((c) => String(c - 1)))}
        value={model.index.pipe(map(String))}
        style={{ '--msp-progress': progress }}
        aria-label="Model"
        onInput={(e) => model.seek(Number(e.currentTarget.value))}
      />
      <Button
        variant="ghost"
        size="sm"
        icon={<Icon name="chevron-left" />}
        title="Previous model"
        onClick={() => step(-1)}
      />
      <PlayButton playing={model.isPlaying} title="Play models" onClick={() => model.togglePlay()} />
      <Button
        variant="ghost"
        size="sm"
        icon={<Icon name="chevron-right" />}
        title="Next model"
        onClick={() => step(1)}
      />
    </div>
  );
}

function PlayButton(props: { playing: ReadonlyState<boolean>; title: string; onClick(): void }) {
  return (
    <Button
      variant="ghost"
      size="sm"
      icon={
        <Show when={props.playing} fallback={() => <Icon name="play" />}>
          {() => <Icon name="pause" />}
        </Show>
      }
      title={props.playing.pipe(map((p) => (p ? 'Pause' : props.title)))}
      pressed={props.playing}
      onClick={() => props.onClick()}
    />
  );
}
