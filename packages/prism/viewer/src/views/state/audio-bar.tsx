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
import { onCleanup } from '@molstar/prism-ui/dom/scope';
import { combineLatest, map } from 'rxjs';
import { type AudioModel, formatTime } from '../../features/snapshots/audio.js';
import { classes } from './_generated/audio-bar.module.css.js';

// A story's audio in the playback card (in place of the browser's own player): play or pause, how far it is as a line
// to click or drag along, the time, stop, and close (the player goes until a command makes another).

export function AudioBar(props: { model: AudioModel }) {
  const { model } = props;
  const fraction = combineLatest([model.time, model.duration]).pipe(map(([t, d]) => (d > 0 ? Math.min(1, t / d) : 0)));
  // Seeking: a press on the line, then a drag along it (the pointer captured); the drag seeks once a frame, not for each
  // move (each seek is the media element's, and a request for streamed audio), and where it was let go.
  let rect: DOMRect | undefined;
  let frame = 0;
  let clientX = 0;
  const seek = () => {
    if (frame) cancelAnimationFrame(frame);
    frame = 0;
    if (rect) model.seek((clientX - rect.left) / (rect.width || 1));
  };
  const release = () => {
    seek();
    rect = undefined;
  };
  onCleanup(() => {
    if (frame) cancelAnimationFrame(frame);
    frame = 0;
  });
  return (
    <div class={classes.bar} role="group" aria-label="Audio">
      <button
        type="button"
        class={classes.play}
        title={model.playing.pipe(map((p) => (p ? 'Pause' : 'Play')))}
        aria-label={model.playing.pipe(map((p) => (p ? 'Pause audio' : 'Play audio')))}
        onClick={() => model.toggle()}
      >
        <Show when={model.playing} fallback={() => <Icon name="play" size={13} />}>
          {() => <Icon name="pause" size={13} />}
        </Show>
      </button>
      <span class={classes.time}>{model.time.pipe(map(formatTime))}</span>
      <div
        class={classes.track}
        role="slider"
        tabIndex={0}
        aria-label="Position"
        aria-valuemin="0"
        aria-valuemax={model.duration.pipe(map((d) => String(Math.round(d))))}
        aria-valuenow={model.time.pipe(map((t) => String(Math.round(t))))}
        aria-valuetext={combineLatest([model.time, model.duration]).pipe(
          map(([t, d]) => `${formatTime(t)} of ${formatTime(d)}`),
        )}
        onPointerDown={(e) => {
          if (e.button !== 0) return;
          const el = e.currentTarget as HTMLElement;
          rect = el.getBoundingClientRect();
          el.setPointerCapture(e.pointerId);
          clientX = e.clientX;
          seek();
        }}
        onPointerMove={(e) => {
          if (!rect) return;
          clientX = e.clientX;
          if (!frame) frame = requestAnimationFrame(seek);
        }}
        onPointerUp={(e) => {
          clientX = e.clientX;
          release();
        }}
        onPointerCancel={() => {
          if (frame) cancelAnimationFrame(frame);
          frame = 0;
          rect = undefined;
        }}
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
        <span class={classes.fill} style={{ transform: fraction.pipe(map((f) => `scaleX(${f})`)) }} />
        <span class={classes.knobTrack} style={{ transform: fraction.pipe(map((f) => `translateX(${f * 100}%)`)) }}>
          <span class={classes.knob} />
        </span>
      </div>
      <span class={classes.time}>{model.duration.pipe(map((d) => (d > 0 ? formatTime(d) : '–:––')))}</span>
      <button type="button" class={classes.action} title="Stop" aria-label="Stop audio" onClick={() => model.stop()}>
        <Icon name="square" size={11} />
      </button>
      <button type="button" class={classes.action} title="Close" aria-label="Close audio" onClick={() => model.close()}>
        <Icon name="x" size={13} />
      </button>
    </div>
  );
}
