/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import {
  colorHex,
  inkForColor,
  InkThreshold,
  luminance,
  parseColor,
  parseColorAlpha,
  themeForLuminance,
} from '../canvas-ink.js';
import { StageModel, type StageTimers } from '../stage-model.js';

// The stage's model in plain Node (stage §1, headless): the ink, the measured inputs with their defaults, the hover
// lines, and the empty mark's wait on fake timers.

describe('canvas ink', () => {
  it('reads luminance as WCAG does', () => {
    expect(luminance(0x000000)).toBe(0);
    expect(luminance(0xffffff)).toBeCloseTo(1, 6);
    expect(luminance(0x808080)).toBeCloseTo(0.2159, 3);
  });

  it('turns at 0.179: dark chrome below, light above', () => {
    expect(themeForLuminance(InkThreshold - 0.001)).toBe('dark');
    expect(themeForLuminance(InkThreshold)).toBe('light');
    expect(inkForColor('#000').theme).toBe('dark');
    expect(inkForColor('#fff').theme).toBe('light');
    // Mid grey (0x76: luminance 0.181) is light, one step darker (0x75: 0.178) dark.
    expect(inkForColor('#767676').theme).toBe('light');
    expect(inkForColor('#757575').theme).toBe('dark');
  });

  it('parses the colors an app sets', () => {
    expect(parseColor('#abc')).toBe(0xaabbcc);
    expect(parseColor('#A1B2C3')).toBe(0xa1b2c3);
    expect(parseColor('#a1b2c3ff')).toBe(0xa1b2c3);
    expect(parseColor('rgb(1, 2, 3)')).toBe(0x010203);
    expect(parseColor('rgba(255 0 0 / 0.5)')).toBe(0xff0000);
    expect(parseColor(' black ')).toBe(0);
    expect(parseColor('linear-gradient(red, blue)')).toBeUndefined();
    expect(colorHex(0x0a0b0c)).toBe('#0a0b0c');
  });

  it('parses names, hsl(), percent channels, alphas and transparent', () => {
    expect(parseColor('navy')).toBe(0x000080);
    expect(parseColor('PapayaWhip')).toBe(0xffefd5);
    expect(parseColor('rgb(50% 0% 0%)')).toBe(0x800000);
    expect(parseColor('rgba(100%, 50%, 0%, 0.5)')).toBe(0xff8000);
    expect(parseColor('hsl(0 0% 5%)')).toBe(0x0d0d0d);
    expect(parseColor('hsla(120, 100%, 25%, 1)')).toBe(0x008000);
    expect(parseColor('hsl(0.5turn 100% 50%)')).toBe(0x00ffff);
    expect(parseColor('#0008')).toBe(0x000000);
    expect(parseColorAlpha('#00000080')?.[1]).toBeCloseTo(0.5, 2);
    expect(parseColorAlpha('rgb(0 0 0 / 0%)')).toEqual([0, 0]);
    expect(parseColorAlpha('transparent')).toEqual([0, 0]);
    expect(parseColor('oklch(0.2 0 0)')).toBeUndefined();
    expect(parseColor('var(--bg)')).toBeUndefined();
    expect(parseColor('constructor')).toBeUndefined();
    expect(parseColor('rgb(1, 2)')).toBeUndefined();
  });

  it('keeps the app’s color; one it cannot read (or a transparent one) is light and unresolved: the view’s', () => {
    expect(inkForColor('#123456')).toEqual({ theme: 'dark', color: '#123456', resolved: true });
    expect(inkForColor('hsl(0 0% 5%)')).toEqual({ theme: 'dark', color: 'hsl(0 0% 5%)', resolved: true });
    expect(inkForColor('navy').theme).toBe('dark');
    expect(inkForColor('oklch(0.2 0 0)')).toEqual({ theme: 'light', color: 'oklch(0.2 0 0)', resolved: false });
    expect(inkForColor('transparent')).toEqual({ theme: 'light', color: 'transparent', resolved: false });
  });
});

describe('StageModel', () => {
  beforeEach(() => {
    vi.useFakeTimers();
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it('derives the ink from the background, a new value only when it changes', () => {
    const stage = new StageModel({ empty: false });
    expect(stage.ink.value).toEqual({ theme: 'light', color: '#ffffff', resolved: true });
    const seen: string[] = [];
    const off = stage.ink.subscribe((ink) => seen.push(ink.theme));
    stage.setBackground('#101418');
    stage.setBackground('#101418');
    expect(stage.ink.value).toEqual({ theme: 'dark', color: '#101418', resolved: true });
    expect(seen.at(-1)).toBe('dark');
    off.unsubscribe();
    stage.dispose();
  });

  it('has headless defaults for what a view measures', () => {
    const stage = new StageModel();
    expect(stage.topInset.value).toBe(0);
    expect(stage.size.value).toEqual({ width: 0, height: 0 });
    // Unknown width: not compact.
    expect(stage.compact.value).toBe(false);
    stage.setSize(320.4, 400);
    expect(stage.size.value).toEqual({ width: 320, height: 400 });
    expect(stage.compact.value).toBe(true);
    stage.setSize(1024, 700);
    expect(stage.compact.value).toBe(false);
    stage.setTopInset(83.6);
    expect(stage.topInset.value).toBe(84);
    // The band hidden: no inset.
    stage.showBand(false);
    expect(stage.band.value).toBe(false);
    expect(stage.topInset.value).toBe(0);
    stage.dispose();
  });

  it('holds the hover lines; an equal list changes nothing', () => {
    const stage = new StageModel();
    const lines = ['ALA 12', 'A'];
    stage.setHover(lines);
    let changes = 0;
    const off = stage.hover.lines.subscribe(() => changes++);
    const before = changes;
    stage.setHover([...lines]);
    expect(changes).toBe(before);
    stage.clearHover();
    expect(stage.hover.lines.value).toEqual([]);
    stage.setHoverPlacement('cursor');
    expect(stage.hoverPlacement.value).toBe('cursor');
    off.unsubscribe();
    stage.dispose();
  });

  it('auto placement (the default): the corner on a roomy stage or unmeasured, the pointer on a narrow or short one', () => {
    const stage = new StageModel();
    expect(stage.hoverPlacement.value).toBe('auto');
    expect(stage.hoverAt.value).toBe('corner');
    stage.setSize(900, 600);
    expect(stage.hoverAt.value).toBe('corner');
    stage.setSize(420, 600);
    expect(stage.hoverAt.value).toBe('cursor');
    stage.setSize(900, 300);
    expect(stage.hoverAt.value).toBe('cursor');
    // Chosen, it stays whatever the size.
    stage.setHoverPlacement('corner');
    expect(stage.hoverAt.value).toBe('corner');
    stage.dispose();
  });

  it('starts empty and quiet; the mark comes after the delay, empty and idle', () => {
    const stage = new StageModel();
    expect(stage.quiet.value).toBe(true);
    expect(stage.mark.value).toBe('gone');
    vi.advanceTimersByTime(699);
    expect(stage.mark.value).toBe('gone');
    vi.advanceTimersByTime(1);
    expect(stage.mark.value).toBe('in');
    expect(stage.quiet.value).toBe(true);
    stage.dispose();
  });

  it('a load never flashes the mark; busy restarts the wait', () => {
    const stage = new StageModel({ empty: false });
    expect(stage.quiet.value).toBe(false);
    stage.setEmpty(true);
    // Cleared to load the next: the canvas still shows until the mark comes.
    expect(stage.quiet.value).toBe(false);
    vi.advanceTimersByTime(500);
    stage.setBusy(true);
    vi.advanceTimersByTime(5000);
    expect(stage.mark.value).toBe('gone');
    stage.setBusy(false);
    vi.advanceTimersByTime(600);
    expect(stage.mark.value).toBe('gone');
    vi.advanceTimersByTime(100);
    expect(stage.mark.value).toBe('in');
    expect(stage.quiet.value).toBe(true);
    // Up: work alone does not take it down.
    stage.setBusy(true);
    vi.advanceTimersByTime(1000);
    expect(stage.mark.value).toBe('in');
    stage.dispose();
  });

  it('fades the mark out only if it showed, then it is gone; the canvas shows', () => {
    const stage = new StageModel();
    vi.advanceTimersByTime(700);
    stage.setEmpty(false);
    expect(stage.mark.value).toBe('out');
    expect(stage.quiet.value).toBe(false);
    vi.advanceTimersByTime(250);
    expect(stage.mark.value).toBe('gone');

    // Never shown: straight to gone.
    const quick = new StageModel();
    vi.advanceTimersByTime(300);
    quick.setEmpty(false);
    expect(quick.mark.value).toBe('gone');
    expect(quick.quiet.value).toBe(false);
    stage.dispose();
    quick.dispose();
  });

  it('runs on injected timers; off, it never goes quiet', () => {
    const waits: number[] = [];
    const stage = new StageModel({
      emptyDelayMs: 50,
      timers: { setTimeout: (_fn, ms) => waits.push(ms), clearTimeout: () => {} },
    });
    expect(waits).toEqual([50]);
    const off = new StageModel({ mark: false });
    expect(off.quiet.value).toBe(false);
    vi.advanceTimersByTime(2000);
    expect(off.mark.value).toBe('gone');
    stage.dispose();
    off.dispose();
  });

  it('a fade out always ends in gone, even when emptied again while busy (a splash assembles again)', () => {
    // Injected timers, run by hand: what the model waits on, in order.
    const pending = new Map<number, { fn: () => void; ms: number }>();
    let next = 0;
    const timers: StageTimers = {
      setTimeout: (fn, ms) => {
        pending.set(++next, { fn, ms });
        return next;
      },
      clearTimeout: (id) => {
        pending.delete(id as number);
      },
    };
    const fire = (ms: number) => {
      const [id, t] = [...pending].find(([, t]) => t.ms === ms)!;
      pending.delete(id);
      t.fn();
    };
    const stage = new StageModel({ timers });
    const phases: string[] = [];
    const off = stage.mark.subscribe((p) => phases.push(p));
    fire(700);
    stage.setEmpty(false);
    expect(stage.mark.value).toBe('out');
    // Emptied again during the fade, while busy (a new load): the fade's timer keeps running.
    stage.setBusy(true);
    stage.setEmpty(true);
    expect([...pending.values()].map((t) => t.ms)).toEqual([250]);
    fire(250);
    expect(stage.mark.value).toBe('gone');
    // Idle: the wait for the mark again, which comes up anew.
    stage.setBusy(false);
    fire(700);
    expect(stage.mark.value).toBe('in');
    // (Its current phase first, if the state tells a new subscriber.)
    expect(phases.slice(phases[0] === 'gone' ? 1 : 0)).toEqual(['in', 'out', 'gone', 'in']);
    // Emptied again during a fade, idle: gone first, then the wait.
    stage.setEmpty(false);
    stage.setEmpty(true);
    fire(250);
    expect(stage.mark.value).toBe('gone');
    expect([...pending.values()].map((t) => t.ms)).toEqual([700]);
    off.unsubscribe();
    stage.dispose();
  });

  it('disposed, waits no more', () => {
    const stage = new StageModel();
    stage.dispose();
    vi.advanceTimersByTime(2000);
    expect(stage.mark.value).toBe('gone');
    expect(vi.getTimerCount()).toBe(0);
  });
});
