/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterEach, beforeEach, describe, expect, it, vi, type MockInstance } from 'vitest';
import { isDebugMode, setDebugMode } from '@molstar/core/util/debug';
import { createGraphicsStats } from '../stats.js';
import { createGraphicsTimer, type TimerResult, type TimerSource } from '../timer.js';

/** Segments take 10 ns each, unless listed in `invalid` */
function createFakeSource(invalid = new Set<number>()) {
  let nextId = 0;
  let current = -1;
  const finished: number[] = [];
  const log: string[] = [];
  const source: TimerSource = {
    begin: () => {
      current = nextId++;
      log.push(`begin ${current}`);
      return current;
    },
    end: () => {
      log.push(`end ${current}`);
      finished.push(current);
      current = -1;
    },
    poll: () => finished.splice(0).map((id) => [id, invalid.has(id) ? undefined : 10] as const),
    clear: () => {
      finished.length = 0;
    },
  };
  return { source, log };
}

function setup(invalid?: Set<number>) {
  const stats = createGraphicsStats();
  const { source, log } = createFakeSource(invalid);
  const timer = createGraphicsTimer(source, stats);
  const draw = () => {
    stats.calls.drawInstanced += 1;
  };
  return { timer, draw, log, stats };
}

function find(results: readonly TimerResult[], label: string): TimerResult | undefined {
  for (const r of results) {
    if (r.label === label) return r;
    const c = find(r.children, label);
    if (c) return c;
  }
}

describe('graphics timer', () => {
  let debugMode = false;
  let warn: MockInstance;
  beforeEach(() => {
    debugMode = isDebugMode;
    setDebugMode(true);
    warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
  });
  afterEach(() => {
    setDebugMode(debugMode);
    warn.mockRestore();
  });

  it('without a source does nothing', () => {
    const stats = createGraphicsStats();
    const timer = createGraphicsTimer(undefined, stats);
    timer.group('A');
    timer.passBegin('drawing-buffer');
    timer.groupEnd('A');
    expect(timer.resolve()).toEqual([]);
  });

  it('sums the passes that begin in a group, nested groups first, then the passes', () => {
    const { timer, draw, log } = setup();
    timer.group('A');
    timer.passBegin('drawing-buffer', 'X');
    draw();
    timer.group('B');
    timer.passBegin('drawing-buffer', 'Y');
    draw();
    timer.passBegin('drawing-buffer', 'Z');
    draw();
    timer.groupEnd('B');
    timer.passBegin('drawing-buffer', 'W');
    draw();
    timer.passEnd();
    timer.groupEnd('A');

    expect(log).toEqual(['begin 0', 'end 0', 'begin 1', 'end 1', 'begin 2', 'end 2', 'begin 3', 'end 3']);
    const results = timer.resolve();
    expect(results.length).toBe(1);
    const a = results[0];
    expect(a.gpuElapsed).toBe(40);
    expect(a.children.map((c) => c.label)).toEqual(['B', 'pass X', 'pass W']);
    const b = find(results, 'B')!;
    expect(b.gpuElapsed).toBe(20);
    expect(b.children.map((c) => [c.label, c.gpuElapsed])).toEqual([
      ['pass Y', 10],
      ['pass Z', 10],
    ]);
    expect(warn).not.toHaveBeenCalled();
  });

  it('times compute passes like render passes', () => {
    const { timer, draw } = setup();
    timer.group('A');
    timer.passBegin('compute', 'scan');
    timer.passBegin('drawing-buffer', 'X');
    draw();
    timer.passBegin('compute');
    timer.passEnd();
    timer.groupEnd('A');

    const a = timer.resolve()[0];
    expect(a.gpuElapsed).toBe(30);
    expect(a.children.map((c) => [c.label, c.gpuElapsed])).toEqual([
      ['pass scan', 10],
      ['pass X', 10],
      ['pass compute', 10],
    ]);
    expect(warn).not.toHaveBeenCalled();
  });

  it('ends the timing of a pass with its innermost group', () => {
    const { timer, draw, log } = setup();
    timer.group('A');
    timer.passBegin('drawing-buffer', 'X');
    draw();
    timer.groupEnd('A');
    // the pass is still open on the GPU, the root resolves
    expect(log).toEqual(['begin 0', 'end 0']);
    expect(timer.resolve()[0].gpuElapsed).toBe(10);
    timer.passEnd();
    expect(log).toEqual(['begin 0', 'end 0']);
  });

  it('gives groups without passes no GPU time', () => {
    const { timer } = setup();
    timer.group('A');
    timer.groupEnd('A');
    const [a] = timer.resolve();
    expect(a.gpuElapsed).toBeNaN();
    expect(a.cpuElapsed).toBeGreaterThanOrEqual(0);
  });

  it('does not time passes begun outside of groups', () => {
    const { timer, draw, log } = setup();
    timer.passBegin('drawing-buffer', 'X');
    draw();
    timer.passEnd();
    expect(log).toEqual([]);
  });

  it('labels passes after their target', () => {
    const { timer } = setup();
    timer.group('A');
    timer.passBegin('drawing-buffer');
    timer.passEnd();
    timer.groupEnd('A');
    expect(timer.resolve()[0].children.map((c) => c.label)).toEqual(['pass drawing-buffer']);
  });

  describe('placement check', () => {
    it('warns about a group inside a pass', () => {
      const { timer, draw } = setup();
      timer.group('A');
      timer.passBegin('drawing-buffer', 'X');
      draw();
      timer.group('B');
      draw();
      timer.groupEnd('B');
      timer.passEnd();
      timer.groupEnd('A');
      expect(warn).toHaveBeenCalledWith("Timer group 'B' begins or ends inside pass 'X'");
      // B has no pass of its own
      expect(find(timer.resolve(), 'B')!.gpuElapsed).toBeNaN();
    });

    it('warns about draws after the group of a pass ended', () => {
      const { timer, draw } = setup();
      timer.group('A');
      timer.group('B');
      timer.passBegin('drawing-buffer', 'X');
      draw();
      timer.groupEnd('B');
      draw();
      timer.passEnd();
      timer.groupEnd('A');
      expect(warn).toHaveBeenCalledWith("Timer group 'A' begins or ends inside pass 'X'");
    });

    it('warns about a group drawing into a pass begun outside of it', () => {
      const { timer, draw } = setup();
      timer.passBegin('drawing-buffer', 'X');
      timer.group('A');
      draw();
      timer.groupEnd('A');
      expect(warn).toHaveBeenCalledWith("Timer group 'A' begins or ends inside pass 'X'");
    });

    it('accepts a group that begins in a pass without draws before it', () => {
      const { timer, draw } = setup();
      timer.passBegin('drawing-buffer', 'X');
      draw();
      // the open pass gets no more draws
      timer.group('A');
      timer.passBegin('drawing-buffer', 'Y');
      draw();
      timer.passEnd();
      timer.groupEnd('A');
      expect(warn).not.toHaveBeenCalled();
    });

    it('is off without debug mode', () => {
      setDebugMode(false);
      const { timer, draw } = setup();
      timer.group('A');
      timer.passBegin('drawing-buffer', 'X');
      timer.group('B');
      draw();
      timer.groupEnd('B');
      timer.passEnd();
      timer.groupEnd('A');
      expect(warn).not.toHaveBeenCalled();
    });
  });

  it('drops results with invalid times', () => {
    const { timer, draw } = setup(new Set([0]));
    timer.group('A');
    timer.passBegin('drawing-buffer', 'X');
    draw();
    timer.passEnd();
    timer.groupEnd('A');
    expect(timer.resolve()).toEqual([]);
  });

  it('captures the draw call statistics', () => {
    const { timer, draw, stats } = setup();
    stats.calls.drawInstanced = 5;
    timer.group('A', { captureStats: true });
    timer.passBegin('drawing-buffer', 'X');
    draw();
    draw();
    timer.passEnd();
    timer.groupEnd('A');
    expect(timer.resolve()[0].calls?.drawInstanced).toBe(2);
    expect(warn).not.toHaveBeenCalled();
  });

  it('throws on unbalanced and duplicate groups', () => {
    const { timer } = setup();
    timer.group('A');
    expect(() => timer.group('A')).toThrow();
    timer.group('B');
    expect(() => timer.groupEnd('A')).toThrow();
    expect(() => timer.groupEnd('C')).toThrow();
  });
});
