/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { lostReleases, stubPointerCapture } from '../../../testing/gestures.js';
import { LineGraphRow } from '../extra-controls.js';
import type { FieldBinding } from '../field-binding.js';

beforeAll(stubPointerCapture);
afterEach(() => {
  document.body.innerHTML = '';
});

describe('LineGraphRow', () => {
  it('a point’s drag (gestures attach once): no listener left by lost releases, and the next drag moves it', () => {
    const start: [number, number][] = [
      [0, 0],
      [0.5, 0.5],
      [1, 1],
    ];
    const points = observableState<[number, number][]>(start);
    const binding: FieldBinding = { value: points.asReadonly, get: () => points.value, set: (v) => points.set(v) };
    const container = document.createElement('div');
    document.body.append(container);
    render(
      () => <LineGraphRow head="Curve" label="Curve" param={PD.LineGraph(start as never)} binding={binding} />,
      container,
    );
    container.querySelector<HTMLButtonElement>('[aria-label=Curve]')!.click();
    const svg = container.querySelector<SVGSVGElement>('svg[viewBox="0 0 240 96"]')!;
    // jsdom lays nothing out: the graph 240 × 96 px at the origin.
    svg.getBoundingClientRect = () => ({ left: 0, top: 0, width: 240, height: 96 }) as DOMRect;
    const r = lostReleases({
      target: () => svg.querySelector('circle[data-index="1"]')!,
      from: { x: 120, y: 48 },
      to: { x: 60, y: 24 },
      reset: () => points.set(start),
      worked: () => points.value.length === 3 && points.value[1][0] === 0.25 && points.value[1][1] === 0.75,
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });
});
