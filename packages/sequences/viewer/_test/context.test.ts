/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Construction and disposal (seqspec:viewer#the-context), the tracked scheduler and the context's guarded fan-out in
// its order (sync before emit). Commands and behaviors are CM1's and IN2's.

import { describe, expect, it } from 'vitest';
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  manualScheduler,
  type SeqIssue,
  SeqViewerConfig,
} from '../src/index.js';
import { trackedScheduler } from '../src/scheduler.js';
import { smallDoc } from './graph-fixtures.js';

describe('construction and disposal', () => {
  it('config: spec items over defaults; set read at next use', () => {
    const ctx = createSeqViewer({ ...DefaultSeqViewerSpec(), config: [[SeqViewerConfig.Viewport.SettleMs, 90]] });
    expect(ctx.config.get(SeqViewerConfig.Viewport.SettleMs)).toBe(90);
    expect(ctx.config.get(SeqViewerConfig.Input.DragSlopPx)).toBe(3);
    ctx.config.set(SeqViewerConfig.Input.DragSlopPx, 5);
    expect(ctx.config.get(SeqViewerConfig.Input.DragSlopPx)).toBe(5);
  });

  it("config by item key: another copy of the viewer's items set the same values (VS-19)", () => {
    // Another bundle's `SeqViewerConfig`: equal items, other objects.
    const other = { ...SeqViewerConfig.Viewport.SettleMs };
    const ctx = createSeqViewer({ ...DefaultSeqViewerSpec(), config: [[other, 90]] });
    expect(ctx.config.get(SeqViewerConfig.Viewport.SettleMs)).toBe(90);
    ctx.config.set({ ...SeqViewerConfig.Input.DragSlopPx }, 7);
    expect(ctx.config.get(SeqViewerConfig.Input.DragSlopPx)).toBe(7);
  });

  it('the manual scheduler fires due timers in order, then frames', () => {
    const s = manualScheduler();
    const log: string[] = [];
    s.setTimer(() => log.push('b'), 20);
    s.setTimer(() => {
      log.push('a');
      s.setTimer(() => log.push('a+5'), 5);
    }, 10);
    s.requestFrame((t) => log.push(`frame ${t}`));
    const c = s.setTimer(() => log.push('c'), 15);
    s.clearTimer(c);
    s.advance(30);
    expect(log).toEqual(['a', 'a+5', 'b', 'frame 30']);
    expect(s.now()).toBe(30);
  });

  it('dispose completes every state and event, clears every timer and frame, and is idempotent', () => {
    const scheduler = manualScheduler();
    const ctx = createSeqViewer(undefined, { scheduler });
    const completed: string[] = [];
    const watch = (name: string, s: { subscribe(o: { complete(): void }): unknown }) =>
      s.subscribe({ complete: () => completed.push(name) });
    for (const [name, s] of Object.entries({
      'registry.version': ctx.registry.version,
      skin: ctx.skin,
      assembled: ctx.graph.assembled,
      shown: ctx.document.shown,
      views: ctx.views.current,
    }))
      watch(name, s);
    const issue = ctx.events.issue as unknown as { subscriberCount: number };
    ctx.events.issue.subscribe(() => {});
    ctx.scheduler.setTimer(() => completed.push('timer fired'), 100);
    ctx.scheduler.requestFrame(() => completed.push('frame fired'));
    expect(scheduler.pendingTimers + scheduler.pendingFrames).toBe(2);
    ctx.dispose();
    ctx.dispose();
    expect(completed.sort()).toEqual(['assembled', 'registry.version', 'shown', 'skin', 'views']);
    expect(scheduler.pendingTimers + scheduler.pendingFrames).toBe(0);
    scheduler.advance(1000);
    expect(completed).not.toContain('timer fired');
    expect(issue.subscriberCount).toBe(0);
    expect(ctx.disposed).toBe(true);
    expect(() => ctx.setSkin(null)).toThrow(/after dispose/);
  });

  it('a timer loop that reschedules after dispose is refused: nothing stays pending', () => {
    const scheduler = manualScheduler();
    const ctx = createSeqViewer(undefined, { scheduler });
    let fired = 0;
    const dwell = () => {
      fired++;
      ctx.scheduler.setTimer(dwell, 10);
    };
    ctx.scheduler.setTimer(dwell, 10);
    scheduler.advance(25);
    expect(fired).toBe(2);
    ctx.dispose();
    dwell(); // a step that was awaiting when the context closed
    expect(scheduler.pendingTimers).toBe(0);
    expect(ctx.scheduler.requestFrame(() => {})).toBe(-1);
    scheduler.advance(100);
    expect(fired).toBe(3);
  });
});

describe('the tracked scheduler', () => {
  it('tracks one callback requested twice, and timers of one callback by due time', () => {
    const port = manualScheduler();
    const s = trackedScheduler(port);
    let [frames, timers] = [0, 0];
    const frame = () => void frames++;
    const timer = () => void timers++;
    s.requestFrame(frame);
    s.requestFrame(frame);
    s.setTimer(timer, 50);
    s.setTimer(timer, 10);
    expect(s.outstanding).toBe(4);
    port.advance(20); // the 10 ms timer, then both frames
    expect([frames, timers, s.outstanding]).toEqual([2, 1, 1]);
    s.clearAll(); // clears the 50 ms one, not the one that fired
    port.advance(100);
    expect([timers, port.pendingTimers, s.closed]).toEqual([1, 0, true]);
  });
});

/** The context's fan-out (internal: the managers' phases add to it). */
interface Fanout {
  sync: ((doc: unknown) => void)[];
  emit: (() => void)[];
  assembled: ((a: unknown) => void)[];
}
const fanoutOf = (ctx: object) => (ctx as { fanout: Fanout }).fanout;

describe("the context's own subscriptions", () => {
  it('sync every manager, then emit; a throw is seqv.state.error and the rest still run', async () => {
    const ctx = createSeqViewer();
    const issues: SeqIssue[] = [];
    ctx.events.issue.subscribe((i) => i.code.startsWith('seqv.') && issues.push(i));
    const order: string[] = [];
    const thrown = new Set<string>();
    const once = (what: string) => () => {
      if (thrown.has(what)) return;
      thrown.add(what);
      throw new Error(what);
    };
    const f = fanoutOf(ctx);
    f.sync.push(
      (d) => order.push(`viewports.sync ${d !== null}`),
      once('layout'),
      () => order.push('ui.keep'),
    );
    f.emit.push(() => order.push('sections'));
    f.assembled.push(once('selection'), (a) => order.push(`focus ${a !== null}`));
    await ctx.setDocument(smallDoc());
    expect(order.slice(0, 3)).toEqual(['viewports.sync true', 'ui.keep', 'sections']);
    expect(order).toContain('focus true');
    expect(issues.map((i) => i.code)).toEqual(['seqv.state.error', 'seqv.state.error']);
    expect(issues.map((i) => i.message)).toEqual([
      'a viewer subscriber threw: layout',
      'a viewer subscriber threw: selection',
    ]);
  });
});
