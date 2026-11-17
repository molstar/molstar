/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import { type Child, Show } from '../dom/index.js';
import { Splash } from '../components/splash/splash.js';
import { classes } from './_generated/empty-stage.module.css.js';
import type { StageModel } from './stage-model.js';

export interface EmptyStageProps {
  model: StageModel;
  /**
   * What shows on the empty stage (stage §4), the app's: rendered as the mark comes up, given whether the stage is
   * compact (a smaller mark). Default: the Mol* glass splash; `() => null` for none.
   */
  mark?: (ctx: { readonly compact: boolean }) => Child;
}

/** The default empty mark: the Mol* wordmark as glass, smaller on a compact stage. */
export const DefaultEmptyMark = (ctx: { readonly compact: boolean }): Child => <Splash scale={ctx.compact ? 2 : 3} />;

/**
 * The empty stage (stage §4): the app's mark (default the Mol* glass) in the middle of the stage while its model's mark
 * is up (empty and idle
 * a while: `StageModel.mark`, whose wait it follows, so a load never flashes it). It fades in, and out faster when
 * something comes, only if it showed; meanwhile the model keeps the canvas quiet. Never in the way of the canvas: only
 * the word takes the pointer.
 */
export function EmptyStage(props: EmptyStageProps) {
  const { model } = props;
  const mark = props.mark ?? DefaultEmptyMark;
  return (
    <Show when={derive(model.mark, (p) => p !== 'gone')}>
      {() => (
        <div class={classes.empty} data-phase={model.mark} data-msui-stage-empty="">
          {mark({ compact: model.compact.value })}
        </div>
      )}
    </Show>
  );
}
