/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import { type Child, Show } from '../dom/index.js';
import { classes } from './_generated/empty-stage.module.css.js';
import type { StageModel } from './stage-model.js';

/** What an empty mark is given: whether the stage is compact (a smaller mark). */
export interface EmptyMarkContext {
  readonly compact: boolean;
}

/** The app's mark on an empty stage (stage §4): rendered as it comes up, and again when `compact` changes. */
export type EmptyMark = (ctx: EmptyMarkContext) => Child;

export interface EmptyStageProps {
  model: StageModel;
  /**
   * What shows on the empty stage (stage §4), the app's: rendered as the mark comes up and again when the stage turns
   * compact or roomy. None by default (the stage is a base: nothing shows, the model's wait and `quiet` still run); the
   * Mol* glass is `GlassEmptyMark` (`./glass-empty-mark.js`).
   */
  mark?: EmptyMark;
}

/**
 * The empty stage (stage §4): the app's mark in the middle of the stage while its model's mark is up (empty and idle a
 * while: `StageModel.mark`, whose wait it follows, so a load never flashes it). It fades in, and out faster when
 * something comes (`StageModel.markFadeOutMs`, as `--msui-stage-mark-out`), only if it showed; meanwhile the model
 * keeps the canvas quiet. Never in the way of the canvas: only the mark takes the pointer. No mark: nothing.
 */
export function EmptyStage(props: EmptyStageProps) {
  const { model, mark } = props;
  if (!mark) return null;
  // One context per size: a mark made anew when the stage turns compact or roomy (a resize changes its size).
  const ctx = derive(model.compact, (compact): EmptyMarkContext => ({ compact }));
  return (
    <Show when={derive(model.mark, (p) => p !== 'gone')}>
      {() => (
        <div
          class={classes.empty}
          data-phase={model.mark}
          data-msui-stage-empty=""
          style={{ '--msui-stage-mark-out': `${model.markFadeOutMs}ms` }}
        >
          <Show when={ctx}>{(c) => mark(c)}</Show>
        </div>
      )}
    </Show>
  );
}
