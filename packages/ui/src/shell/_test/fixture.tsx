/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { vi } from 'vitest';
import { domView, flush, getScope, onCleanup, render, type ViewScope } from '../../dom/index.js';
import { ShellModel, type ShellOptions, type ShellViewDef } from '../shell-model.js';
import { ShellView, type ShellViewProps } from '../shell-view.js';

// A shell mounted in jsdom with small inline views (no story files): each view records its scope (to check it is
// paused while not shown), its cleanups, and holds an input whose typed value tells whether its body was kept.

export interface Probe {
  scope?: ViewScope;
  mounts: number;
  cleanups: number;
}

export function probeView(id: string, extra: Partial<ShellViewDef> = {}): { def: ShellViewDef; probe: Probe } {
  const probe: Probe = { mounts: 0, cleanups: 0 };
  const def: ShellViewDef = {
    id,
    title: id[0].toUpperCase() + id.slice(1),
    icon: 'file',
    render: domView(() => {
      probe.mounts++;
      probe.scope = getScope();
      onCleanup(() => probe.cleanups++);
      return (
        <div data-probe={id}>
          <input data-input={id} aria-label={`${id} draft`} />
          <button type="button" data-button={id}>
            {id} action
          </button>
        </div>
      );
    }),
    ...extra,
  };
  return { def, probe };
}

export interface Mounted {
  readonly shell: ShellModel;
  readonly root: HTMLElement;
  readonly probes: Record<string, Probe>;
  dispose(): void;
  q<T extends Element = HTMLElement>(selector: string): T | null;
  qa<T extends Element = HTMLElement>(selector: string): T[];
  rail(): HTMLElement | null;
  railButton(title: string): HTMLButtonElement | null;
  /** The docked node of a group (its wrapper in the tree: `data-shown`, `data-folded`). */
  groupNode(id: string): HTMLElement | null;
  frame(ref: string): HTMLElement | null;
  win(ref: string): HTMLElement | null;
  main(): HTMLElement;
}

export function mountShell(
  views: { def: ShellViewDef; probe: Probe }[],
  options: ShellOptions = {},
  props: Partial<ShellViewProps> = {},
): Mounted {
  const shell = new ShellModel(options);
  shell.registerViews(views.map((v) => v.def));
  const container = document.createElement('div');
  container.className = 'msui-root';
  document.body.append(container);
  const dispose = render(
    // The ⋯ button on (the tests reach the layout menu by it); `viewMenu: false` for the default.
    () => <ShellView shell={shell} main={() => <div data-app-main="">scene</div>} viewMenu {...props} />,
    container,
  );
  flush();
  const probes = Object.fromEntries(views.map((v) => [v.def.id, v.probe]));
  const q = <T extends Element = HTMLElement>(s: string) => container.querySelector<T>(s);
  const qa = <T extends Element = HTMLElement>(s: string) => [...container.querySelectorAll<T>(s)];
  return {
    shell,
    root: container,
    probes,
    dispose: () => {
      dispose();
      container.remove();
    },
    q,
    qa,
    rail: () => q('nav'),
    railButton: (title) => q<HTMLButtonElement>(`nav button[aria-label="${title}"]`),
    groupNode: (id) => q(`[data-group="${id}"]`)?.closest<HTMLElement>('[data-node]') ?? null,
    frame: (ref) => q(`[data-msui-view="${ref}"]`),
    win: (ref) => q(`[data-msui-view="${ref}"]`)?.closest<HTMLElement>('[data-window]') ?? null,
    main: () => q<HTMLElement>('[data-msui-shell-main]')!,
  };
}

export function key(target: EventTarget, k: string, init: KeyboardEventInit = {}): KeyboardEvent {
  const e = new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true, ...init });
  target.dispatchEvent(e);
  flush();
  return e;
}

/** A press (pointerdown, then a click on the same element): what dismissal listens to. */
export function press(el: Element): void {
  el.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, cancelable: true, composed: true }));
  flush();
}

export const microtasks = async () => {
  await Promise.resolve();
  await Promise.resolve();
  flush();
};

/**
 * jsdom has no ResizeObserver: a stub that records what is observed, and `resize(el, size)` that sets the element's
 * client size and calls its observers (as a layout would).
 */
export function stubResizeObserver(): { resize(el: Element, size: { width?: number; height?: number }): void } {
  const observed = new Map<Element, Set<ResizeObserverCallback>>();
  class Stub {
    readonly #fn: ResizeObserverCallback;
    readonly #els = new Set<Element>();
    constructor(fn: ResizeObserverCallback) {
      this.#fn = fn;
    }
    observe(el: Element) {
      this.#els.add(el);
      let s = observed.get(el);
      if (!s) observed.set(el, (s = new Set()));
      s.add(this.#fn);
    }
    unobserve(el: Element) {
      observed.get(el)?.delete(this.#fn);
      this.#els.delete(el);
    }
    disconnect() {
      for (const el of this.#els) observed.get(el)?.delete(this.#fn);
      this.#els.clear();
    }
  }
  vi.stubGlobal('ResizeObserver', Stub);
  return {
    resize(el, size) {
      if (size.width !== undefined) Object.defineProperty(el, 'clientWidth', { value: size.width, configurable: true });
      if (size.height !== undefined)
        Object.defineProperty(el, 'clientHeight', { value: size.height, configurable: true });
      for (const fn of [...(observed.get(el) ?? [])]) fn([], {} as ResizeObserver);
      flush();
    },
  };
}
