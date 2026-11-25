/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's scroll past its cap (_spec/README.md "Input", The one rail; core's "Render", Section scroll): the viewer's
// own, composed wherever the cap is (no behavior needed: a minimal viewer scrolls its sections too). Its top in px of
// the scrolled content: a tracks section's viewport `top`, a text's anchor moved by px through its placed lines
// (`Placement.scrollTopAt`). Plain functions over the last frame's section (a wheel allocates nothing) and the section rail's gesture
// (a drag on its window, a press beside it centring there, Esc puts back the top the press found; pan only, session
// state: no history step).

import type { SectionLayout, ViewLayout, WindowPart } from '@molstar/seqx-core/render';
import { cancelFrameCallback, requestFrameCallback } from '@molstar/ui/dom/frame';
import type { ViewModel } from '@molstar/seqx-viewer/models';

/** The section of the last frame by its address that scrolls past its cap (its `sectionRail`), else undefined. */
export function scrolledSection(layout: ViewLayout | undefined, section: string): SectionLayout | undefined {
  const ss = layout?.sections;
  if (ss) for (let i = 0; i < ss.length; i++) if (ss[i].section.address === section && ss[i].sectionRail) return ss[i];
  return undefined;
}

/** The px a section's scrolled content can move: its rail's total less its window. */
export function sectionMax(s: SectionLayout): number {
  const r = s.sectionRail;
  return r ? Math.max(0, r.total - (r.end - r.start)) : 0;
}

/** A section's top now (px), from the stored viewport (moves before a frame add up), clamped to its content. */
export function sectionTop(view: ViewModel, s: SectionLayout): number {
  const v = view.get(s.section.address);
  const p = s.placement;
  const t =
    s.layout.name === 'text' ? (v && p.scrollTopAt ? p.scrollTopAt(v.start) : (p.top ?? 0)) : (v?.top ?? p.top ?? 0);
  return Math.max(0, Math.min(sectionMax(s), t));
}

/** A section's top at `top` (px), clamped to `[0, max]`; false when it was there. */
export function setSectionTop(view: ViewModel, s: SectionLayout, top: number): boolean {
  const now = sectionTop(view, s);
  const next = Math.max(0, Math.min(sectionMax(s), top));
  if (Math.abs(next - now) < 1e-6) {
    // a stored top past the content (it shrank) comes back within it
    const stored = view.get(s.section.address)?.top;
    if (s.layout.name !== 'text' && stored !== undefined && stored !== next) view.setTop(s.section.address, next);
    return false;
  }
  if (s.layout.name === 'text') return view.panPx(s.section.address, 0, next - now);
  view.setTop(s.section.address, next);
  return true;
}

/**
 * The wheel's vertical delta over a section past its cap, coalesced to the next frame as the navigation's wheel is
 * (`navigate.coalesceWheel`): a wheel event allocates nothing (the model's viewport record is published once a frame).
 * One pending section at a time; a wheel over another applies the first's at once.
 */
export class SectionWheel {
  private _section = '';
  private _dy = 0;
  private _frame: number | undefined;
  private readonly _run = (): void => void this.apply();

  constructor(
    private readonly view: ViewModel,
    private readonly layout: () => ViewLayout | undefined,
  ) {}

  /** A wheel's `dy` px over `section`: false at its end with nothing pending that way (the wheel goes to the page). */
  wheel(section: string, dy: number): boolean {
    const s = dy ? scrolledSection(this.layout(), section) : undefined;
    if (!s) return false;
    if (this._dy && section !== this._section) this.apply();
    const now = sectionTop(this.view, s) + this._dy;
    const next = Math.max(0, Math.min(sectionMax(s), now + dy));
    if (Math.abs(next - now) < 1e-6) return false;
    this._section = section;
    this._dy += next - now;
    this._frame ??= requestFrameCallback(this._run);
    return true;
  }

  /** The pending scroll applied now (the frame's; a key or a gesture before it reads the model). */
  apply(): boolean {
    if (this._frame !== undefined) cancelFrameCallback(this._frame);
    this._frame = undefined;
    const dy = this._dy;
    this._dy = 0;
    const s = dy ? scrolledSection(this.layout(), this._section) : undefined;
    return s ? setSectionTop(this.view, s, sectionTop(this.view, s) + dy) : false;
  }

  dispose(): void {
    if (this._frame !== undefined) cancelFrameCallback(this._frame);
    this._frame = undefined;
    this._dy = 0;
  }
}

/**
 * Page y `y0` to `y1` (the last frame's: a band's row) shown in a section past its cap by the least scroll (the
 * keyboard cursor's reveal); nothing where the section does not scroll.
 */
export function revealInSection(
  view: ViewModel,
  layout: ViewLayout | undefined,
  section: string,
  y0: number,
  y1: number,
): void {
  const s = scrolledSection(layout, section);
  if (!s) return;
  const p = s.placement;
  const from = p.scrollY ?? p.box.y;
  const to = p.scrollEnd || p.box.y + p.box.height;
  const now = sectionTop(view, s);
  if (y0 < from) setSectionTop(view, s, now - (from - y0));
  else if (y1 > to) setSectionTop(view, s, now + (y1 - to));
}

/**
 * After a layout: a tracks section's stored `top` kept within its content while it scrolls (its tracks hidden or shrunk
 * while scrolled; the layout clamps what it draws by too). Allocates nothing.
 */
export function clampSectionTops(view: ViewModel, layout: ViewLayout): void {
  const ss = layout.sections;
  for (let i = 0; i < ss.length; i++) {
    const s = ss[i];
    // only where it scrolls: a top on a section that shows all of it is read by nothing (a story step's)
    if (s.layout.name === 'text' || !s.sectionRail) continue;
    const a = s.section.address;
    const stored = view.get(a)?.top;
    if (stored === undefined || stored === 0) continue;
    const max = sectionMax(s);
    if (stored > max) view.setTop(a, max);
  }
}

/** The section rail's gesture: a drag on its window scrolls by the px the pointer travelled along it. */
export class SectionRailGesture {
  private _press: { section: string; at0: number; start: number; part: WindowPart; saved: number } | undefined;

  constructor(
    private readonly view: ViewModel,
    private readonly layout: () => ViewLayout | undefined,
  ) {}

  get pressed(): boolean {
    return this._press !== undefined;
  }

  /** A press at `at` (px along the scrolled content); false when the section does not scroll. */
  press(section: string, part: WindowPart, at: number): boolean {
    const s = scrolledSection(this.layout(), section);
    if (!s || !Number.isFinite(at)) return false;
    const top = sectionTop(this.view, s);
    this._press = { section, at0: at, start: top, part, saved: top };
    return true;
  }

  dragTo(at: number): boolean {
    const p = this._press;
    const s = p && p.part !== 'outside' ? scrolledSection(this.layout(), p.section) : undefined;
    if (!p || !s || !Number.isFinite(at)) return false;
    return setSectionTop(this.view, s, p.start + at - p.at0);
  }

  release(): void {
    this._press = undefined;
  }

  /** Esc: the top the press found (by `setTop`: a follower's too). */
  cancel(): void {
    const p = this._press;
    this._press = undefined;
    const s = p && scrolledSection(this.layout(), p.section);
    if (p && s) setSectionTop(this.view, s, p.saved);
  }

  /** A click beside the window (a press that never moved): the window centred there; on it nothing. */
  click(section: string, part: WindowPart, at: number, detail = 1): boolean {
    this._press = undefined;
    const s = scrolledSection(this.layout(), section);
    if (!s || detail > 1 || part !== 'outside' || !Number.isFinite(at)) return false;
    const r = s.sectionRail!;
    return setSectionTop(this.view, s, at - (r.end - r.start) / 2);
  }
}
