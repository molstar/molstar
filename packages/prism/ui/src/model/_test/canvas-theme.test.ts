/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Subject } from 'rxjs';
import { observableState } from '../../reactive/observable-state.js';
import {
  type CanvasBackgroundProps,
  canvasInk,
  canvasTheme,
  canvasThemeOf,
  parseColor,
  themeForColor,
} from '../canvas-theme.js';

const solid = (backgroundColor: number): CanvasBackgroundProps => ({ renderer: { backgroundColor } });

describe('canvas theme', () => {
  it('is dark on a dark background and light on a light one', () => {
    expect(themeForColor(0x000000)).toBe('dark');
    expect(themeForColor(0x14181f)).toBe('dark');
    expect(themeForColor(0xffffff)).toBe('light');
    // Mid gray: dark text still has the contrast.
    expect(themeForColor(0x777777)).toBe('light');
    expect(themeForColor(0x555555)).toBe('dark');
  });

  it('averages a gradient over the background color by its opacity', () => {
    const gradient = (opacity: number): CanvasBackgroundProps => ({
      renderer: { backgroundColor: 0x000000 },
      postprocessing: {
        background: {
          variant: { name: 'radialGradient', params: { centerColor: 0xffffff, edgeColor: 0xdddddd, opacity } },
        },
      },
    });
    expect(canvasThemeOf(gradient(1))).toBe('light');
    expect(canvasThemeOf(gradient(0.1))).toBe('dark');
  });

  it('uses the background color under an image, and the fallback for a transparent canvas', () => {
    expect(canvasThemeOf({ ...solid(0x000000), postprocessing: { background: { variant: { name: 'image' } } } })).toBe(
      'dark',
    );
    expect(canvasThemeOf({ ...solid(0x000000), transparentBackground: true })).toBe('light');
    expect(canvasThemeOf({ ...solid(0xffffff), transparentBackground: true }, 'dark')).toBe('dark');
  });

  it('follows the props, emitting only changes', () => {
    const props = new Subject<CanvasBackgroundProps>();
    const seen: string[] = [];
    canvasTheme(props).subscribe((t) => seen.push(t));
    props.next(solid(0xffffff));
    props.next(solid(0xeeeeee));
    props.next(solid(0x000000));
    props.next(solid(0x111111));
    props.next(solid(0xffffff));
    expect(seen).toEqual(['light', 'dark', 'light']);
  });

  it('parses hex and rgb colors', () => {
    expect(parseColor('#fff')).toBe(0xffffff);
    expect(parseColor('14181f')).toBe(0x14181f);
    expect(parseColor('rgba(20, 24, 31, 0.5)')).toBe(0x14181f);
    expect(parseColor('nope')).toBeUndefined();
  });

  it('keeps an ink state in step with the props, distinct on change', () => {
    const props = observableState<CanvasBackgroundProps>(solid(0xffffff));
    const ink = canvasInk(props.asReadonly);
    const seen: string[] = [];
    const sub = ink.subscribe((i) => seen.push(`${i.theme}${i.color}`));
    expect(ink.value).toEqual({ theme: 'light', color: '#ffffff' });
    props.set(solid(0x14181f));
    props.set({ ...solid(0x14181f), transparentBackground: false });
    expect(ink.value).toEqual({ theme: 'dark', color: '#14181f' });
    expect(seen).toEqual(['light#ffffff', 'dark#14181f']);
    props.set({ ...solid(0x000000), transparentBackground: true });
    expect(ink.value).toEqual({ theme: 'light', color: '#ffffff' });
    sub.unsubscribe();
    ink.dispose();
  });
});
