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
  tap,
  timer,
} from 'rxjs';
import { MolstarWordmark } from '../../components/icon/wordmark.js';
import { Show } from '../../dom/control-flow.js';
import { classes } from './_generated/empty-canvas.module.css.js';
import { type PrismFieldOptions, PrismSplash, type SplashKind } from './prism-splash.js';

/**
 * How long the canvas stays empty, the viewer idle, before the mark (and its splash) shows: clearing a scene to load the
 * next never flashes it, and any work restarts the wait (a long load or a session opening shows nothing meanwhile).
 */
const Delay = 700;
/** How often the random splash is the wave (`quantum-field.ts`); the rest, the light's particles. */
const QuantumShare = 0.8;
/** How long it takes to fade out when something loads (faster than it came). */
const FadeOut = 250;

type Phase = 'in' | 'out' | 'gone';

/**
 * The Mol* wordmark in the middle of an empty canvas, once it has been empty and the viewer idle a while (work restarts
 * the wait): it fades in as the logo's colors wipe in along the outline (and stay); when something loads it fades out,
 * faster (and if it never showed, nothing shows). With `splash`, light through the Mol* star comes with it
 * (`PrismSplash`), as particles or as a wave, picked at random each time it shows (`field`). Once shown, work alone does not hide it (a load that fills the canvas does). Never in the way of the
 * canvas.
 */
export function EmptyCanvas(props: {
  empty: Observable<boolean>;
  /** The splash with the mark (on by default); its options, to tune it (the light's: they pick the light). */
  splash?: boolean | Partial<PrismFieldOptions>;
  /**
   * Which splash: the light (`prism-field.ts`), the light as a wave (`quantum-field.ts`), or either, picked at random
   * each time it shows (the default: the wave four times in five, `QuantumShare`).
   */
  field?: SplashKind | 'random';
  /** The viewer is working (a load, a task): the wait for the mark starts over while it is. */
  busy?: Observable<boolean>;
  /**
   * Something loaded while the mark showed (the canvas had been empty and the viewer idle): the host reveals it
   * gradually. Never for a load the mark never showed for (a scene cleared to load the next, work in progress).
   */
  onReveal?: () => void;
}) {
  const phase = combineLatest([props.empty.pipe(distinctUntilChanged()), props.busy ?? of(false)]).pipe(
    switchMap(
      ([empty, busy]): Observable<Phase> =>
        !empty
          ? concat(of('out' as const), timer(FadeOut).pipe(map(() => 'gone' as const)))
          : // Working: nothing changes (a shown mark stays; one not yet shown waits), and the wait starts over after.
            busy
            ? NEVER
            : timer(Delay).pipe(map(() => 'in' as const)),
    ),
    startWith<Phase>('gone'),
    // Fading out only what showed.
    scan((prev: Phase, next: Phase): Phase => (next === 'out' && prev !== 'in' ? 'gone' : next)),
    distinctUntilChanged(),
    tap((p) => p === 'out' && props.onReveal?.()),
    shareReplay({ bufferSize: 1, refCount: true }),
  );
  const splash = props.splash ?? true;
  const field = (): SplashKind => {
    // Options are the light's: with them, the light.
    if (typeof splash === 'object') return 'light';
    const f = props.field ?? 'random';
    return f === 'random' ? (Math.random() < QuantumShare ? 'quantum' : 'light') : f;
  };
  return (
    <Show when={phase.pipe(map((p) => p !== 'gone'))}>
      {() => {
        // The word alone: with the splash, its star is the prism.
        const mark = (
          <MolstarWordmark scale={3} flash="reveal" live={!splash} prism={!splash} class={classes.mark} />
        ) as Element;
        return (
          <div class={classes.empty} data-phase={phase}>
            {splash ? (
              <div class={classes.splash}>
                <PrismSplash
                  field={field()}
                  mark={() => mark}
                  // The viewport: the pointer over the canvas (which keeps it) moves the light too.
                  pointer={() => mark.closest('[data-msp-viewport]') ?? mark.parentElement}
                  options={splash === true ? {} : splash}
                />
              </div>
            ) : null}
            {mark}
          </div>
        );
      }}
    </Show>
  );
}
