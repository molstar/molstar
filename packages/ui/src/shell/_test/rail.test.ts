/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { fitRail, RailPad, RailSectionGap, RailSteps, railEntries, railStep } from '../rail.js';
import type { RailItem, ShellViewDef } from '../shell-model.js';

const view = (id: string, extra: Partial<ShellViewDef> = {}) => ({
  ref: id,
  def: { id, title: id, icon: 'file', render: { mount: () => () => {} }, ...extra } as ShellViewDef,
});
const item = (id: string, extra: Partial<RailItem> = {}): RailItem => ({
  id,
  title: id,
  icon: 'file',
  run() {},
  ...extra,
});

describe('railEntries', () => {
  it('orders head, strip, item groups, foot, then the items without a group', () => {
    const es = railEntries(
      [view('help', { place: 'foot' }), view('a'), view('home', { place: 'head' }), view('b', { place: 'strip' })],
      [item('theme'), item('cmd', { group: 'g' })],
    );
    expect(es.map((e) => e.id)).toEqual(['view:home', 'view:a', 'view:b', 'item:cmd', 'view:help', 'item:theme']);
    expect(es.map((e) => e.section)).toEqual(['head', 'strip', 'strip', 'group:g', 'foot', 'tools']);
  });

  it('a corner item comes last (under the Hide button), whatever its group, and stays when the rail runs short', () => {
    const es = railEntries(
      [view('home', { place: 'head' }), view('a'), view('b'), view('help', { place: 'foot' })],
      [item('commands', { corner: true, group: 'g' }), item('theme')],
    );
    expect(es.map((e) => e.id)).toEqual(['view:home', 'view:a', 'view:b', 'view:help', 'item:theme', 'item:commands']);
    expect(es.at(-1)?.section).toBe('corner');
    // Too short for all of them: the others leave first.
    const fit = fitRail(es, RailSteps[RailSteps.length - 1].button * 3 + 2 * RailPad + 3 * RailSectionGap);
    expect(fit.hidden.size).toBeGreaterThan(0);
    expect(fit.hidden.has('item:commands')).toBe(false);
  });
});

describe('fitRail', () => {
  const entries = railEntries(
    [
      view('home', { place: 'head' }),
      view('a'),
      view('b', { priority: 100 }),
      view('c'),
      view('help', { place: 'foot' }),
    ],
    [],
  );

  it('leaves nothing out before the rail is measured or when everything fits', () => {
    expect(fitRail(entries, 0).hidden.size).toBe(0);
    expect(fitRail(entries, 1000).hidden.size).toBe(0);
    expect(fitRail(entries, 1000).step.name).toBe('normal');
  });

  it('puts the highest priority under More first, then the later ones (the foot before the strip)', () => {
    const small = RailSteps.find((s) => s.name === 'small')!;
    // Room for the padding, the rail's own button, More, three sections' gaps and three entries.
    const h = 2 * RailPad + small.button * 5 + 2 * RailSectionGap;
    const fit = fitRail(entries, h);
    expect(fit.step.name).toBe('small');
    // b (priority 100) first, then help (the last, by position), then c.
    expect([...fit.hidden]).toEqual(['view:b', 'view:help']);
    const tighter = fitRail(entries, h - small.button);
    expect([...tighter.hidden]).toEqual(['view:b', 'view:help', 'view:c']);
  });

  it('steps by height at the steps’ minimums', () => {
    for (const s of RailSteps) expect(railStep(s.min).name).toBe(s.name);
    expect(railStep(RailSteps[0].min - 1).name).toBe(RailSteps[1].name);
  });
});

describe('RailSteps and shell.module.css', () => {
  const css = readFileSync(resolve(process.cwd(), 'packages/ui/src/shell/shell.module.css'), 'utf8');
  const unit = 4;
  const px = (expr: string) => {
    const m = /calc\((\d+) \* var\(--msui-unit\)\)/.exec(expr);
    if (!m) throw new Error(`not a unit expression: ${expr}`);
    return Number(m[1]) * unit;
  };

  it('the normal step’s button is --msui-shell-button on .shell', () => {
    const m = /\.shell\s*\{[^}]*--msui-shell-button:\s*([^;]+);/.exec(css);
    expect(m).not.toBeNull();
    expect(px(m![1])).toBe(RailSteps[0].button);
  });

  it('every smaller step is a @container shell (max-height) rule one px under the step above’s minimum', () => {
    const rules = [
      ...css.matchAll(/@container shell \(max-height: (\d+)px\)\s*\{[^{]*\{\s*--msui-shell-button:\s*([^;]+);/g),
    ].map((m) => ({ max: Number(m[1]), button: px(m[2]) }));
    expect(rules).toHaveLength(RailSteps.length - 1);
    for (let i = 1; i < RailSteps.length; i++) {
      expect(rules[i - 1]).toEqual({ max: RailSteps[i - 1].min - 1, button: RailSteps[i].button });
    }
  });

  it('a captioned button is the step’s button plus the caption', () => {
    const m = /--msui-shell-caption:\s*([^;]+);/.exec(css);
    const caption = px(m![1]);
    for (const s of RailSteps) expect(s.captioned).toBe(s.button + caption);
  });
});
