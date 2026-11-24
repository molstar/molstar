/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The square's drag frame allocates nothing once warm: over the ring with its cells (`rondo:square`), a drag of its
// run's box in progress, a step a pointer move and the square's frame (its preview, no patch until the release, its
// cells not laid out again). By hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import type { SquareNavigator } from '@molstar/seqx-core/render/navigator-square';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { gripPoint, squareOf } from '../square-gestures.js';
import { ringStoryView } from '../views.js';

describe('the square over the ring', () => {
  it('a drag frame of its box (a move and the frame)', async ({ task }) => {
    await loadGenome();
    const v = await openExample(ringStoryView('ring'), {
      delayMs: 0,
      wait: 3000,
      fixtures: true,
      extensions: [rondoViewer],
    });
    await v.execute('section.include', { address: '/views/main/ring', include: ['chr4', 'chr5', 'chr6'] });
    v.draw();
    const n = squareOf(v);
    const [x, y] = gripPoint(n, 'body')!;
    n.down(x, y);
    n.move(x + 20, y + 20);
    n.draw();
    // a drag step: the gesture's snap of the pointer (set in the square's `pointer`, as the model's move sets it) and its
    // preview, then the square's frame of it (its base kept, the window pass drawn again), its input reused
    const priv = n as unknown as { renderer: SquareNavigator; setGesture: { drag(s: SquareNavigator): void } };
    const square = priv.renderer;
    const gesture = priv.setGesture;
    const move = quickCase(task, 'square drag move', (i) => {
      const d = ((i * 7) % 60) - 30;
      square.pointer[0] = x + d;
      square.pointer[1] = y + d;
      gesture.drag(square);
    }).r;
    const input = {
      size: { width: n.width, height: n.height },
      dpr: 2,
      settled: true,
      view: v.frame!.layout,
      moves: {},
      pointerX: x,
      pointerY: y,
    };
    square.frame(input as never);
    const frame = quickCase(task, 'square drag frame', (i) => {
      const d = ((i * 7) % 60) - 30;
      square.pointer[0] = x + d;
      square.pointer[1] = y + d;
      gesture.drag(square);
      square.frame(input as never);
    }).r;
    // the model's own step, as the canvas drives it: `NavigatorModel.move` in the drag, then its `draw` (its frame's
    // input, the square's frame, the readout)
    // (the points made once, held as objects, as a pointer event's numbers are: the step measures the model, not the
    // boxing of a number the test computes)
    const points: unknown[] = Array.from({ length: 60 }, (_, d) => [x + d - 30 + 0.5, y + d - 30 + 0.5]).flat();
    points.push('held');
    const model = quickCase(task, 'square drag, the model (move and draw)', (i) => {
      const k = 2 * ((i * 7) % 60);
      n.move(points[k] as number, points[k + 1] as number);
      n.draw();
    }).r;
    console.info(
      `square drag move: ${move.bytesPerStep.toFixed(0)} B; frame: ${frame.bytesPerStep.toFixed(0)} B; model: ${model.bytesPerStep.toFixed(0)} B`,
    );
    n.cancel();
    v.retire();
    expect(move.bytesPerStep, move.note).toBeLessThanOrEqual(inject('perfThreshold'));
    expect(frame.bytesPerStep, frame.note).toBeLessThanOrEqual(inject('perfThreshold'));
    expect(model.bytesPerStep, model.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
