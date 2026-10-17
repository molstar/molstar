/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { lostReleases, stubPointerCapture } from '@molstar/prism-ui/testing/gestures';
import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import type { MeasurementPoint, MeasurementsModel } from '../../../model/measurements.js';
import type { SceneNode } from '../../../model/scene.js';
import { MeasurementsAdd } from '../measurements-add.js';

// The measurement points: dragged by a grip to reorder (`drag`, dom/scoped.ts: attached once, every end path ends it).

beforeAll(stubPointerCapture);
afterEach(() => {
  document.body.innerHTML = '';
});

describe('the measurement points', () => {
  it('drags a point by its grip: no listener left by lost releases, and the next drag moves it', () => {
    const moved: [string, number][] = [];
    const model = {
      points: observableState<readonly MeasurementPoint[]>(['a', 'b', 'c'].map((id) => ({ id, label: `Point ${id}` })))
        .asReadonly,
      movePointTo: (id: string, index: number) => moved.push([id, index]),
      kinds: [],
      added: observableState<readonly SceneNode[]>([]).asReadonly,
    } as unknown as MeasurementsModel;
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <MeasurementsAdd model={model} />, container);
    // jsdom lays nothing out: rows 20 px high, one under another.
    const items = [...container.querySelectorAll<HTMLElement>('[data-id]')];
    expect(items).toHaveLength(3);
    items.forEach((el, i) => (el.getBoundingClientRect = () => ({ top: 20 * i, height: 20 }) as DOMRect));
    const r = lostReleases({
      target: () => items[0].querySelector<HTMLElement>('[data-tip="Drag to reorder"]')!,
      from: { x: 5, y: 5 },
      to: { x: 5, y: 55 },
      reset: () => (moved.length = 0),
      worked: () =>
        moved.length === 1 &&
        moved[0][0] === 'a' &&
        moved[0][1] === 2 &&
        !items.some((el) => el.hasAttribute('data-dragging') || el.hasAttribute('data-drop')),
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });
});
