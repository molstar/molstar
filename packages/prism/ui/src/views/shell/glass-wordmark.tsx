/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { WordmarkGeometry } from '../../components/icon/wordmark.js';
import { onActive, onCleanup } from '../../dom/scope.js';
import { classes } from './_generated/glass-wordmark.module.css.js';
import { Glass } from './glass.js';

/** How far the pieces travel, as a factor of the empty canvas's: a heading has little room. */
const Travel = 0.35;
/** Room around the word for the pieces (px at 1×). */
const Margin = 22;

function prefersReducedMotion(): boolean {
  return typeof matchMedia === 'function' && matchMedia('(prefers-reduced-motion: reduce)').matches;
}

/**
 * The Mol* wordmark as glass (`glass.ts`), without the light: Home's heading. It assembles from its pieces each time it
 * shows, its pieces lean under the pointer, and a press shatters it for a moment before it snaps back. It draws only
 * while something moves: whole and still, nothing is repainted until the pointer comes. Under reduced motion, or where
 * there is nothing to draw it with, the word whole.
 */
export function GlassWordmark(props: { scale?: number; class?: string }): Node {
  const scale = Math.max(1, Math.round(props.scale ?? 1));
  const width = WordmarkGeometry.width * scale;
  const height = WordmarkGeometry.height * scale;
  const margin = Margin * scale;
  const canvas = (
    <canvas
      class={classes.canvas}
      aria-hidden="true"
      style={{
        left: `${-margin}px`,
        top: `${-margin}px`,
        width: `${width + 2 * margin}px`,
        height: `${height + 2 * margin}px`,
      }}
    />
  ) as HTMLCanvasElement;
  const root = (
    <span
      class={props.class ? `${classes.wordmark} ${props.class}` : classes.wordmark}
      style={{ width: `${width}px`, height: `${height}px` }}
      role="img"
      aria-label="Mol*"
    >
      {canvas}
    </span>
  ) as HTMLElement;

  const ctx = canvas.getContext('2d');
  const still = prefersReducedMotion() || typeof Path2D !== 'function' || !ctx;
  const glass = still ? undefined : new Glass(23, Travel);
  // The word's drawing origin in the canvas (its viewBox runs from -1): px.
  const origin = margin + scale;
  let dpr = 1;
  let gradient: CanvasGradient | string = 'currentColor';
  let edge = 'currentColor';

  /** The logo's gradient and the edges' color, from the theme where it shows; the canvas sized for the screen. */
  const prepare = () => {
    if (!ctx) return;
    dpr = Math.min(2, globalThis.devicePixelRatio || 1);
    canvas.width = Math.round((width + 2 * margin) * dpr);
    canvas.height = Math.round((height + 2 * margin) * dpr);
    const color = (name: string) => {
      canvas.style.color = `var(${name})`;
      return getComputedStyle(canvas).color;
    };
    const g = ctx.createLinearGradient(0, 0, WordmarkGeometry.width - 2, 0);
    g.addColorStop(0, color('--msp-color-brand-1'));
    g.addColorStop(0.35, color('--msp-color-brand-2'));
    g.addColorStop(0.6, color('--msp-color-brand-3'));
    g.addColorStop(1, color('--msp-color-brand-4'));
    gradient = g;
    edge = color('--msp-color-text-subtle');
    canvas.style.color = '';
  };

  const draw = () => {
    if (!ctx) return;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    if (glass) {
      glass.draw(ctx, origin, origin, scale, gradient, edge, 1);
      return;
    }
    // Still: the word whole, as the glass draws it at rest.
    ctx.translate(origin, origin);
    ctx.scale(scale, scale);
    ctx.strokeStyle = gradient;
    ctx.lineWidth = WordmarkGeometry.stroke(scale);
    ctx.lineJoin = 'miter';
    for (let i = 0; i < WordmarkGeometry.paths.length; i++) ctx.stroke(new Path2D(WordmarkGeometry.paths[i]));
  };

  // Frames only while something moves.
  let frame = 0;
  let last = 0;
  let active = true;
  const tick = (now: number) => {
    frame = 0;
    if (!glass || !active) return;
    const dt = last ? Math.min(0.05, (now - last) / 1000) : 1 / 60;
    last = now;
    glass.step(dt);
    draw();
    if (!glass.resting) frame = requestAnimationFrame(tick);
    else last = 0;
  };
  const wake = () => {
    if (!glass || !active || frame) return;
    frame = requestAnimationFrame(tick);
  };

  // The pointer, in drawing units: the box is read when it comes in (not per move: no layout reads at pointer rate).
  let box: DOMRect | undefined;
  const at = (e: PointerEvent) => {
    box ??= root.getBoundingClientRect();
    return { u: (e.clientX - box.left) / scale - 1, v: (e.clientY - box.top) / scale - 1 };
  };
  if (glass) {
    root.addEventListener('pointerenter', () => {
      box = root.getBoundingClientRect();
    });
    root.addEventListener('pointermove', (e) => {
      const p = at(e);
      glass.setPointer(p.u, p.v);
      wake();
    });
    root.addEventListener('pointerleave', () => {
      glass.clearPointer();
      box = undefined;
      wake();
    });
    root.addEventListener('pointerdown', (e) => {
      const p = at(e);
      if (glass.phase !== 'broken') glass.shatter(p.u, p.v);
      wake();
    });
  }

  // Each time it shows (the Home tab opened): it assembles. Hidden: no frames.
  onActive((a) => {
    active = a;
    if (!a) {
      if (frame) cancelAnimationFrame(frame);
      frame = 0;
      last = 0;
      return;
    }
    // Sized and colored once it is in the page (on its first frame).
    requestAnimationFrame(() => {
      if (!active) return;
      prepare();
      if (glass) {
        glass.scatter();
        wake();
      } else draw();
    });
  });
  onCleanup(() => {
    active = false;
    if (frame) cancelAnimationFrame(frame);
  });
  return root;
}
