/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { Markdown } from '../../components/markdown/markdown.js';
import { Show } from '../../dom/control-flow.js';
import type { SnapshotStepperModel, TrajectoryControlsModel } from '../../model/playback.js';
import { classes } from './_generated/playback.module.css.js';

/**
 * The snapshot panel at the top left of the canvas: the stepper and the current snapshot's description. Nothing is
 * stacked under it, so a description that grows or collapses moves nothing else.
 */
export function SnapshotPanel(props: { model: SnapshotStepperModel }) {
  return (
    <div class={classes.cluster}>
      <SnapshotStepper model={props.model} />
    </div>
  );
}

function SnapshotStepper(props: { model: SnapshotStepperModel }) {
  const { model } = props;
  const visible = model.count.pipe(map((c) => c > 0));
  const counter = combineLatest([model.index, model.count]).pipe(map(([i, c]) => `${i + 1} / ${c}`));
  return (
    <Show when={visible}>
      {() => (
        <section class={classes.card} aria-label="Snapshots">
          <div class={classes.bar}>
            <Button
              variant="ghost"
              size="sm"
              icon={<Icon name="chevron-left" />}
              title="Previous snapshot"
              onClick={() => model.previous()}
            />
            <span class={classes.counter}>{counter}</span>
            <Button
              variant="ghost"
              size="sm"
              icon={<Icon name="chevron-right" />}
              title="Next snapshot"
              onClick={() => model.next()}
            />
            <span class={classes.title} title={model.title}>
              {model.title}
            </span>
            <Button
              variant="ghost"
              size="sm"
              icon={
                <Show when={model.isPlaying} fallback={() => <Icon name="play" />}>
                  {() => <Icon name="pause" />}
                </Show>
              }
              title={model.isPlaying.pipe(map((p) => (p ? 'Pause' : 'Play snapshots')))}
              pressed={model.isPlaying}
              onClick={() => model.togglePlay()}
            />
          </div>
          <Show
            when={model.transition.pipe(
              map((t) => t !== undefined),
              distinctUntilChanged(),
            )}
          >
            {() => (
              <div
                class={classes.progress}
                style={{ transform: model.transition.pipe(map((t) => `scaleX(${t ?? 0})`)) }}
              />
            )}
          </Show>
          <Show when={model.description}>
            {(source) => (
              <div class={classes.description}>
                <Markdown source={source} onLink={(href) => model.openLink(href)} />
              </div>
            )}
          </Show>
        </section>
      )}
    </Show>
  );
}

/**
 * The timeline at the bottom right of the canvas, like a media player: model index of a trajectory with a scrubber
 * (later also animations). Shown for multi-model trajectories only.
 */
export function Timeline(props: { model: TrajectoryControlsModel }) {
  const { model } = props;
  const visible = model.count.pipe(map((c) => c > 1));
  return (
    <Show when={visible}>
      {() => (
        <section class={`${classes.card} ${classes.bar} ${classes.timeline}`} aria-label="Trajectory">
          <Icon name="trajectory" size={16} />
          <span class={classes.counter}>
            {combineLatest([model.index, model.count]).pipe(map(([i, c]) => `${i + 1} / ${c}`))}
          </span>
          <Button
            variant="ghost"
            size="sm"
            icon={<Icon name="skip-back" />}
            title="First model"
            onClick={() => model.first()}
          />
          <Button
            variant="ghost"
            size="sm"
            icon={<Icon name="step-back" />}
            title="Previous model"
            onClick={() => model.previous()}
          />
          <input
            class={classes.scrubber}
            type="range"
            min="0"
            max={model.count.pipe(map((c) => String(c - 1)))}
            value={model.index.pipe(map(String))}
            aria-label="Model"
            onInput={(e) => model.seek(Number(e.currentTarget.value))}
          />
          <Button
            variant="ghost"
            size="sm"
            icon={<Icon name="step-forward" />}
            title="Next model"
            onClick={() => model.next()}
          />
          <Button
            variant="ghost"
            size="sm"
            icon={<Icon name="skip-forward" />}
            title="Last model"
            onClick={() => model.last()}
          />
          <Button
            variant="ghost"
            size="sm"
            icon={
              <Show when={model.isPlaying} fallback={() => <Icon name="play" />}>
                {() => <Icon name="pause" />}
              </Show>
            }
            title="Play models"
            pressed={model.isPlaying}
            onClick={() => model.togglePlay()}
          />
        </section>
      )}
    </Show>
  );
}
