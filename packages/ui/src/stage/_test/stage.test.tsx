/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import axe from 'axe-core';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { type Child, DomHost, flush, render, ViewScope } from '../../dom/index.js';
import { fakeContext, stubPath2D } from '../../components/logo/_test/fake-canvas.js';
import { classes } from '../_generated/stage.module.css.js';
import { BandPlaceholder, structureLines, StagePlaceholder } from '../placeholders.js';
import { Stage } from '../stage.js';
import { StageModel } from '../stage-model.js';
import { StageToolbar } from '../toolbar.js';

// The stage (stage §2): its layers and their pointer rules, the ink, the band's inset measured into the model, the
// empty mark, and axe over its states.

const host = DomHost.requestAnimationFrame;
let queue: ((t: number) => void)[] = [];
let now = 0;
const frames = (n = 1) => {
  for (let i = 0; i < n && queue.length; i++) {
    now += 16;
    const run = queue;
    queue = [];
    for (const f of run) f(now);
  }
};

beforeEach(() => {
  stubPath2D();
  vi.stubGlobal('matchMedia', () => ({ matches: false }));
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockReturnValue(fakeContext().ctx as never);
  DomHost.requestAnimationFrame = (f) => {
    queue.push(f);
    return 0;
  };
});

let current: { dispose(): void; stage: StageModel } | undefined;
afterEach(() => {
  current?.dispose();
  current?.stage.dispose();
  current = undefined;
  document.body.innerHTML = '';
  frames(200);
  DomHost.requestAnimationFrame = host;
  queue = [];
  vi.useRealTimers();
});

function mount(
  stage: StageModel,
  slots: () => {
    top?: Child;
    lead?: Child;
    status?: Child;
    chips?: Child;
    corner?: Child;
    toolbarShown?: boolean | State<boolean>;
  } = () => ({}),
  scope = new ViewScope(),
) {
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(
    () => (
      <Stage
        model={stage}
        canvas={<StagePlaceholder onHover={(i) => stage.setHover(i === undefined ? [] : structureLines(i))} />}
        {...slots()}
      />
    ),
    container,
    scope,
  );
  current = { dispose, stage };
  const root = container.querySelector<HTMLElement>('[data-msui-stage]')!;
  return {
    root,
    canvas: root.querySelector<HTMLElement>('[data-msui-stage-canvas]')!,
    ink: root.querySelector<HTMLElement>(`.${classes.ink}`)!,
    bottom: root.querySelector<HTMLElement>(`.${classes.bottom}`)!,
  };
}

const wheel = (el: Element, ctrlKey: boolean) => {
  const e = new WheelEvent('wheel', { ctrlKey, bubbles: true, cancelable: true });
  el.dispatchEvent(e);
  return e.defaultPrevented;
};

describe('Stage', () => {
  it('paints the canvas color and themes its chrome by the ink', () => {
    const stage = new StageModel({ empty: false });
    const { root, ink } = mount(stage);
    expect(root.style.getPropertyValue('--msui-stage-color')).toBe('#ffffff');
    expect(ink.getAttribute('data-msui-theme')).toBe('light');
    stage.setBackground('#0b0e14');
    flush();
    expect(root.style.getPropertyValue('--msui-stage-color')).toBe('#0b0e14');
    expect(ink.getAttribute('data-msui-theme')).toBe('dark');
  });

  it('layers the canvas, then the chrome; the slots hold the app’s content', () => {
    const stage = new StageModel({ empty: false });
    const { root, canvas, ink } = mount(stage, () => ({
      top: <p>band</p>,
      lead: <p>lead</p>,
      status: <p>status</p>,
      chips: <p>chips</p>,
      corner: <p>corner</p>,
    }));
    expect([...root.children]).toEqual([canvas, ink]);
    expect(ink.querySelector(`.${classes.band}`)?.textContent).toBe('band');
    expect(ink.querySelector(`.${classes.lead}`)?.textContent).toBe('lead');
    expect(ink.querySelector(`.${classes.status}`)?.textContent).toBe('status');
    expect(ink.querySelector(`.${classes.chips}`)?.textContent).toBe('chips');
    expect(ink.querySelector(`.${classes.corner}`)?.textContent).toContain('corner');
  });

  it('shows its toolbar as the app says (a value or a state: a pane’s status), hidden faded out', () => {
    const stage = new StageModel({ empty: false });
    const shown = new State(false);
    const { root } = mount(stage, () => ({ corner: <p>tools</p>, toolbarShown: shown }));
    expect(root.getAttribute('data-toolbar')).toBe('hidden');
    expect(root.querySelector(`.${classes.tools}`)?.textContent).toBe('tools');
    shown.set(true);
    flush();
    expect(root.hasAttribute('data-toolbar')).toBe(false);
    current?.dispose();
    stage.dispose();
    const always = mount(new StageModel({ empty: false }), () => ({ corner: <p>tools</p> }));
    expect(always.root.hasAttribute('data-toolbar')).toBe(false);
  });

  it('holds a pinch over the canvas and the bottom bar only; Safari’s gestures on the root', () => {
    const stage = new StageModel({ empty: false });
    const { root, canvas, bottom } = mount(stage, () => ({
      top: <p class="band">band</p>,
      status: <p class="status">s</p>,
    }));
    expect(wheel(canvas.firstElementChild!, true)).toBe(true);
    expect(wheel(canvas, false)).toBe(false);
    expect(wheel(bottom.querySelector('.status')!, true)).toBe(true);
    expect(wheel(root.querySelector('.band')!, true)).toBe(false);
    const g = new Event('gesturestart', { cancelable: true });
    root.dispatchEvent(g);
    expect(g.defaultPrevented).toBe(true);
  });

  it('measures the band into the model and writes the inset for the chrome below it', () => {
    const stage = new StageModel({ empty: false });
    const { root } = mount(stage, () => ({ top: <BandPlaceholder /> }));
    // The band's content as a layout would report it, its fade not counted.
    resize(root.querySelector(`.${classes.band}`)!, 52);
    expect(stage.topInset.value).toBe(52);
    expect(root.style.getPropertyValue('--msui-stage-top-inset')).toBe('52px');
    // Hidden: paused and inert, no inset.
    stage.showBand(false);
    flush();
    const top = root.querySelector<HTMLElement>(`.${classes.top}`)!;
    expect(top.hidden).toBe(true);
    expect(top.hasAttribute('inert')).toBe(true);
    expect(root.style.getPropertyValue('--msui-stage-top-inset')).toBe('0px');
  });

  it('measures its size into the model (compact below 480)', () => {
    const stage = new StageModel({ empty: false });
    const { root } = mount(stage);
    resize(root, 300, 320);
    expect(stage.size.value).toEqual({ width: 320, height: 300 });
    expect(stage.compact.value).toBe(true);
    expect(root.hasAttribute('data-compact')).toBe(true);
  });

  it('shows the mark after the wait, quiets the canvas, and takes it out when something comes', () => {
    vi.useFakeTimers();
    const stage = new StageModel();
    const { root, canvas } = mount(stage);
    expect(canvas.hasAttribute('data-quiet')).toBe(true);
    expect(root.querySelector('[data-msui-stage-empty]')).toBeNull();
    vi.advanceTimersByTime(700);
    flush();
    const empty = root.querySelector<HTMLElement>('[data-msui-stage-empty]')!;
    expect(empty.getAttribute('data-phase')).toBe('in');
    expect(empty.querySelector('[data-msui-splash]')).not.toBeNull();
    stage.setEmpty(false);
    flush();
    expect(canvas.hasAttribute('data-quiet')).toBe(false);
    expect(empty.getAttribute('data-phase')).toBe('out');
    vi.advanceTimersByTime(250);
    flush();
    expect(root.querySelector('[data-msui-stage-empty]')).toBeNull();
  });

  it('feeds the hover label from the canvas, at the corner above the toolbar', async () => {
    const stage = new StageModel({ empty: false });
    const { root } = mount(stage, () => ({ corner: <p>tools</p> }));
    await Promise.resolve();
    const bead = root.querySelector('[data-index="3"]')!;
    bead.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    expect(stage.hover.lines.value).toEqual(structureLines(3));
    bead.dispatchEvent(new PointerEvent('pointermove', { bubbles: true, clientX: 5, clientY: 5 }));
    flush();
    const label = root.querySelector<HTMLElement>(`.${classes.corner} [role=tooltip]`)!;
    const title = structureLines(3)[0];
    expect(label.textContent).toContain(typeof title === 'string' ? title : title.runs[0].text);
    expect(label.hasAttribute('aria-hidden')).toBe(false);
    stage.setHoverPlacement('cursor');
    flush();
    expect(root.querySelector(`.${classes.corner} [role=tooltip]`)).toBeNull();
    expect(root.querySelector('[role=tooltip][data-placement="cursor"]')).not.toBeNull();
  });
});

/** Calls the observers of `el` with an entry of its content box, as a layout would. */
function resize(el: Element, height: number, width = 800) {
  const entry = { target: el, contentRect: { width, height } } as unknown as ResizeObserverEntry;
  for (const fn of [...(observers.get(el) ?? [])]) fn([entry], {} as ResizeObserver);
  flush();
}

/** jsdom has no ResizeObserver: a stub keeping the observers by element. */
const observers = new Map<Element, Set<ResizeObserverCallback>>();
beforeEach(() => {
  observers.clear();
  class Stub {
    readonly #fn: ResizeObserverCallback;
    constructor(fn: ResizeObserverCallback) {
      this.#fn = fn;
    }
    observe(el: Element) {
      let s = observers.get(el);
      if (!s) observers.set(el, (s = new Set()));
      s.add(this.#fn);
    }
    unobserve() {}
    disconnect() {
      for (const s of observers.values()) s.delete(this.#fn);
    }
  }
  vi.stubGlobal('ResizeObserver', Stub);
});

// axe over the stage's states: WCAG 2.2 A and AA and best practices, less what jsdom cannot evaluate (no layout, no
// painted colors) and the page-level rules (the stage is mounted in a bare container).
const Off = [
  'color-contrast',
  'color-contrast-enhanced',
  'target-size',
  'scrollable-region-focusable',
  'document-title',
  'html-has-lang',
  'html-lang-valid',
  'landmark-one-main',
  'page-has-heading-one',
  'region',
  'bypass',
];

async function violations(root: Element): Promise<string> {
  const rules: axe.RuleObject = {};
  for (const id of Off) rules[id] = { enabled: false };
  const r = await axe.run(root, {
    runOnly: { type: 'tag', values: ['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa', 'wcag22aa', 'best-practice'] },
    rules,
    iframes: false,
    resultTypes: ['violations'],
  });
  return r.violations
    .map((v) => `${v.id}: ${v.help}\n${v.nodes.map((n) => `  ${n.html.slice(0, 140)}`).join('\n')}`)
    .join('\n');
}

describe('the stage is accessible (axe)', () => {
  const tools = () => (
    <StageToolbar
      groups={[
        {
          id: 'camera',
          label: 'Camera',
          items: [
            { id: 'reset', icon: 'refresh-cw', title: 'Reset the camera', run: () => {} },
            {
              id: 'mode',
              icon: 'camera',
              title: 'Projection',
              menu: [
                { id: 'p', label: 'Perspective', run: () => {} },
                { id: 'o', label: 'Orthographic', run: () => {} },
              ],
            },
          ],
        },
      ]}
    />
  );
  const states: [string, () => StageModel, (s: StageModel) => Record<string, Child>][] = [
    ['empty, the mark up', () => new StageModel(), () => ({})],
    [
      'loaded, toolbar and status',
      () => new StageModel({ empty: false }),
      () => ({ corner: tools(), status: <span>Ready</span> }),
    ],
    [
      'a band on top, dark',
      () => new StageModel({ empty: false, background: '#000000' }),
      () => ({ top: <BandPlaceholder selection={[4, 9]} /> }),
    ],
    [
      'hover at the cursor',
      () => new StageModel({ empty: false, hoverPlacement: 'cursor' }),
      () => ({ corner: tools() }),
    ],
  ];
  for (const [name, make, slots] of states) {
    it(name, async () => {
      vi.useFakeTimers();
      const stage = make();
      const { root } = mount(stage, () => slots(stage));
      vi.advanceTimersByTime(700);
      stage.setHover(structureLines(2));
      flush();
      vi.useRealTimers();
      expect(await violations(root)).toBe('');
    });
  }
});
