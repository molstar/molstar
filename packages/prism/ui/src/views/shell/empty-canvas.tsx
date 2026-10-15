/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  concat,
  distinctUntilChanged,
  map,
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

/** How long the canvas stays empty before the mark shows: clearing a scene to load the next never flashes it. */
const Delay = 700;
/** How long it takes to fade out when something loads (faster than it came). */
const FadeOut = 250;

type Phase = 'in' | 'out' | 'gone';

/**
 * The Mol* wordmark in the middle of an empty canvas, a while after it became empty: it fades in, then the logo's colors
 * wipe in along the outline and stay; when something loads it fades out, faster (and if it never showed, nothing shows). Never in the
 * way of the canvas.
 */
export function EmptyCanvas(props: { empty: Observable<boolean> }) {
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
  return (
    <Show when={phase.pipe(map((p) => p !== 'gone'))}>
      {() => (
        <div class={classes.empty} data-phase={phase}>
          <MolstarWordmark scale={3} flash="reveal" live class={classes.mark} />
        </div>
      )}
    </Show>
  );
}
