/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import { effect, frameState } from '@molstar/prism-ui/dom/scope';
import { drag } from '@molstar/prism-ui/dom/scoped';
import { derive } from '@molstar/prism-ui/reactive/state';
import { type AudioModel, formatTime } from '../../features/snapshots/audio.js';
import { classes } from './_generated/audio-bar.module.css.js';

// A story's audio in the playback card (in place of the browser's own player): play or pause, how far it is as a line
// to click or drag along, the time, stop, and close (the player goes until a command makes another).

export function AudioBar(props: { model: AudioModel }) {
  const { model } = props;
  const fraction = derive([model.time, model.duration], (t, d) => (d > 0 ? Math.min(1, t / d) : 0));
  // Seeking: a press on the line, then a drag along it (the pointer captured; `drag`, dom/scoped.ts: its listeners
  // attached once, every way it ends ends it); the drag seeks once a frame, not for each move (each seek is the media
  // element's, and a request for streamed audio), and where it was let go.
  let rect: DOMRect | undefined;
  let clientX = 0;
  const seek = () => {
    if (rect) model.seek((clientX - rect.left) / (rect.width || 1));
  };
  // The moves, counted: a frame state, so a drag seeks once a frame (its value is only a tick).
  const moves = frameState(0);
  let n = 0;
  effect(moves, seek);
  const track = (
    <div
      class={classes.track}
      role="slider"
      tabIndex={0}
      aria-label="Position"
      aria-valuemin="0"
      aria-valuemax={derive(model.duration, (d) => String(Math.round(d)))}
      aria-valuenow={derive(model.time, (t) => String(Math.round(t)))}
      aria-valuetext={derive([model.time, model.duration], (t, d) => `${formatTime(t)} of ${formatTime(d)}`)}
      onKeyDown={(e) => {
        const d = model.duration.value;
        if (!(d > 0)) return;
        const step = e.shiftKey ? 10 : 5;
        if (e.key === 'ArrowRight') model.seek((model.time.value + step) / d);
        else if (e.key === 'ArrowLeft') model.seek((model.time.value - step) / d);
        else if (e.key === ' ') model.toggle();
        else return;
        e.preventDefault();
      }}
    >
      <span class={classes.fill} style={{ transform: derive(fraction, (f) => `scaleX(${f})`) }} />
      <span class={classes.knobTrack} style={{ transform: derive(fraction, (f) => `translateX(${f * 100}%)`) }}>
        <span class={classes.knob} />
      </span>
    </div>
  ) as HTMLElement;
  drag(track, {
    start: (e) => {
      if (e.button !== 0) return undefined;
      rect = track.getBoundingClientRect();
      clientX = e.clientX;
      seek();
      return rect;
    },
    move: (e) => {
      clientX = e.clientX;
      moves.offer(++n);
    },
    end: (_r, how, e) => {
      // Let go: a seek where it was (a step still waiting for its frame is then a no-op); cancelled: none.
      if (how === 'up' && e) {
        clientX = e.clientX;
        seek();
      }
      rect = undefined;
    },
  });
  return (
    <div class={classes.bar} role="group" aria-label="Audio">
      <button
        type="button"
        class={classes.play}
        title={derive(model.playing, (p) => (p ? 'Pause' : 'Play'))}
        aria-label={derive(model.playing, (p) => (p ? 'Pause audio' : 'Play audio'))}
        onClick={() => model.toggle()}
      >
        <Show when={model.playing} fallback={() => <Icon name="play" size={13} />}>
          {() => <Icon name="pause" size={13} />}
        </Show>
      </button>
      <span class={classes.time}>{derive(model.time, formatTime)}</span>
      {track}
      <span class={classes.time}>{derive(model.duration, (d) => (d > 0 ? formatTime(d) : '–:––'))}</span>
      <button type="button" class={classes.action} title="Stop" aria-label="Stop audio" onClick={() => model.stop()}>
        <Icon name="square" size={11} />
      </button>
      <button type="button" class={classes.action} title="Close" aria-label="Close audio" onClick={() => model.close()}>
        <Icon name="x" size={13} />
      </button>
    </div>
  );
}
