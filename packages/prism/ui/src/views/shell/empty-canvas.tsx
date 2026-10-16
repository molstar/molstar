/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  combineLatest,
  concat,
  distinctUntilChanged,
  map,
  NEVER,
  type Observable,
  of,
  scan,
  shareReplay,
  startWith,
  switchMap,
  timer,
} from 'rxjs';
import { MolstarWordmark } from '../../components/icon/wordmark.js';
import { Show } from '../../dom/control-flow.js';
import { classes } from './_generated/empty-canvas.module.css.js';
import { type PrismFieldOptions, PrismSplash } from './prism-splash.js';

/** How long the canvas stays empty before the mark shows: clearing a scene to load the next never flashes it. */
const Delay = 700;
/** How long it takes to fade out when something loads (faster than it came). */
const FadeOut = 250;
/** How long the mark has shown, the viewer idle, before the splash fades in under it. */
const SplashDelay = 2500;

type Phase = 'in' | 'out' | 'gone';

/**
 * The Mol* wordmark in the middle of an empty canvas, a while after it became empty: it fades in as the logo's colors
 * wipe in along the outline (and stay); when something loads it fades out, faster (and if it never showed, nothing shows).
 * With `splash`, once the mark has shown and the viewer has been idle a while longer, light through Prism's mark fades in
 * under it (`PrismSplash`). Never in the way of the canvas.
 */
export function EmptyCanvas(props: {
  empty: Observable<boolean>;
  /** The splash under the mark (on by default); its options, to tune it. */
  splash?: boolean | Partial<PrismFieldOptions>;
  /** The viewer is working (a load, a task): the splash waits for it to be idle. */
  busy?: Observable<boolean>;
}) {
  const phase = props.empty.pipe(
    distinctUntilChanged(),
    switchMap(
      (empty): Observable<Phase> =>
        empty
          ? timer(Delay).pipe(map(() => 'in' as const))
          : concat(of('out' as const), timer(FadeOut).pipe(map(() => 'gone' as const))),
    ),
    startWith<Phase>('gone'),
    // Fading out only what showed.
    scan((prev: Phase, next: Phase): Phase => (next === 'out' && prev !== 'in' ? 'gone' : next)),
    distinctUntilChanged(),
    shareReplay({ bufferSize: 1, refCount: true }),
  );
  const splash = props.splash ?? true;
  // Idle on the empty canvas: the mark in, nothing working, for a while.
  const idle = combineLatest([phase, props.busy ?? of(false)]).pipe(
    // Fading out (something loaded), it stays as it was and goes with the mark.
    switchMap(([p, busy]) =>
      p === 'out' ? NEVER : p === 'in' && !busy ? timer(SplashDelay).pipe(map(() => true)) : of(false),
    ),
    distinctUntilChanged(),
  );
  return (
    <Show when={phase.pipe(map((p) => p !== 'gone'))}>
      {() => {
        const mark = (<MolstarWordmark scale={3} flash="reveal" live class={classes.mark} />) as Element;
        return (
          <div class={classes.empty} data-phase={phase}>
            {splash ? (
              <Show when={idle}>
                {() => (
                  <div class={classes.splash}>
                    <PrismSplash mark={() => mark} options={splash === true ? {} : splash} />
                  </div>
                )}
              </Show>
            ) : null}
            {mark}
          </div>
        );
      }}
    </Show>
  );
}
