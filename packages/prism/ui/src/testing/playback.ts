/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  PeriodicBoxMode,
  SnapshotStepperModel,
  SnapshotTransitionModel,
  TrajectoryControlsModel,
} from '../model/playback.js';
import { observableState } from '../reactive/observable-state.js';

const Story = [
  {
    title: 'Overview',
    description: 'Cytochrome **P450 3A4** with its heme cofactor.\n\nSee the [active site](#active-site).',
  },
  {
    title: 'Active site',
    description:
      '### Active site\n- Heme iron coordinated by `Cys442`\n- Substrate channel above the heme\n\nBack to the [overview](#overview).',
  },
  { title: 'Surface', description: undefined as string | undefined },
];

/**
 * Snapshots: the three-step 1TQN story, or `count` generated ones (a long session: the position becomes a list).
 * `post` reports links that are not snapshots.
 */
export function mockSnapshots(post: (text: string) => void, count?: number): SnapshotStepperModel {
  const entries =
    count === undefined
      ? Story
      : Array.from({ length: count }, (_, i) => ({
          title:
            ['Overview', 'Active site', 'Heme pocket', 'Substrate channel', 'Surface', 'Ligand contacts'][i % 6] +
            (i >= 6 ? ` ${Math.floor(i / 6) + 1}` : ''),
          description: i % 3 === 2 ? undefined : `Step ${i + 1} of the tour.`,
        }));
  const index = observableState(0);
  const isPlaying = observableState(false);
  const title = observableState<string | undefined>(entries[0].title);
  const description = observableState<string | undefined>(entries[0].description);
  let timer: ReturnType<typeof setInterval> | undefined;
  const go = (i: number) => {
    const n = (i + entries.length) % entries.length;
    index.set(n);
    title.set(entries[n].title);
    description.set(entries[n].description);
  };
  return {
    count: observableState(entries.length).asReadonly,
    index: index.asReadonly,
    title: title.asReadonly,
    titles: observableState(entries.map((e) => e.title)).asReadonly,
    description: description.asReadonly,
    isPlaying: isPlaying.asReadonly,
    previous: () => go(index.value - 1),
    next: () => go(index.value + 1),
    goTo: go,
    togglePlay: () => {
      isPlaying.set(!isPlaying.value);
      clearInterval(timer);
      if (isPlaying.value) timer = setInterval(() => go(index.value + 1), 2500);
    },
    openLink: (href) => {
      const i = entries.findIndex((e) => `#${e.title.toLowerCase().replace(/ /g, '-')}` === href);
      if (i >= 0) go(i);
      else post(`Link: ${href}`);
    },
  };
}

/** A multi-model trajectory of `count` models. */
export function mockTrajectory(count = 20, label = '1TQN'): TrajectoryControlsModel {
  const index = observableState(0);
  const isPlaying = observableState(false);
  let timer: ReturnType<typeof setInterval> | undefined;
  return {
    label: observableState(label).asReadonly,
    count: observableState(count).asReadonly,
    index: index.asReadonly,
    isPlaying: isPlaying.asReadonly,
    seek: (i) => index.set(Math.max(0, Math.min(count - 1, i))),
    togglePlay: () => {
      isPlaying.set(!isPlaying.value);
      clearInterval(timer);
      if (isPlaying.value) timer = setInterval(() => index.set((index.value + 1) % count), 150);
    },
  };
}

/** A simulation's frames in a periodic box: the trajectory with the box menu. */
export function mockPeriodicTrajectory(count = 100, label = 'AdK in water'): TrajectoryControlsModel {
  const mode = observableState<PeriodicBoxMode | undefined>('whole');
  const images = observableState(false);
  return {
    ...mockTrajectory(count, label),
    unit: observableState<'models' | 'frames'>('frames').asReadonly,
    periodic: {
      mode: mode.asReadonly,
      set: (m) => mode.set(m),
      images: { shown: images.asReadonly, set: (on) => images.set(on) },
    },
  };
}

/** A snapshot's transition of `duration` ms, played in frames of `step` ms. */
export function mockTransition(duration = 2000, step = 100): SnapshotTransitionModel {
  const time = observableState(0);
  const isPlaying = observableState(false);
  let timer: ReturnType<typeof setInterval> | undefined;
  const stop = () => {
    clearInterval(timer);
    isPlaying.set(false);
  };
  return {
    duration: observableState(duration).asReadonly,
    step: observableState(step).asReadonly,
    time: time.asReadonly,
    isPlaying: isPlaying.asReadonly,
    seek: (t) => {
      stop();
      time.set(Math.max(0, Math.min(duration, t)));
    },
    togglePlay: () => {
      if (isPlaying.value) return stop();
      time.set(0);
      isPlaying.set(true);
      timer = setInterval(() => {
        const next = time.value + step;
        time.set(Math.min(duration, next));
        if (next >= duration) stop();
      }, step);
    },
  };
}
