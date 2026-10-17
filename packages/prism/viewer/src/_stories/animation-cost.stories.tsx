/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { For } from '@molstar/prism-ui/dom/control-flow';
import { render } from '@molstar/prism-ui/dom/render';
import { onCleanup, RuntimeProbe } from '@molstar/prism-ui/dom/scope';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { DockId } from '@molstar/prism-ui';
import { story } from '@molstar/prism-ui/testing/index';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { createPrismViewer } from '../page.js';
import { ExampleFiles } from '../testing/real-viewer.js';
import type { PrismViewer } from '../viewer.js';

const meta: Meta = { title: 'Internals/Performance/Animation' };
export default meta;

// What the UI costs while a structure animates (styling.md §5, rendering.md §1): 1GRM, an NMR entry of several models,
// played at a constant rate (`trajectory.play`, a fixed time a pass), with the panels in one of five setups. Per frame:
// binding writes and `watch` renders (the runtime's probe), DOM mutations by area, and the long animation frames'
// script time split into the UI's (Prism, Seq*) and Mol*'s. Mounted outside Storybook's root (its highlight addon
// re-scans styles on every change there), a11y off. A driver adds a CPU profile (`__mspAnim`; watch-bench).

/** Models a second: the playback's constant rate (a pass takes at least a second: 1GRM's 5 models, one second). */
const ModelsPerSecond = 5;

export type AnimationSetup = 'closed' | 'structure' | 'state-tree' | 'sequence' | 'menu';
const Setups: readonly AnimationSetup[] = ['closed', 'structure', 'state-tree', 'sequence', 'menu'];
const Docks: readonly DockId[] = ['left', 'right', 'top', 'bottom'];

/** The areas mutations are counted by: the first whose element holds the mutated node. */
const Areas: readonly [string, string][] = [
  ['scene tree', '[role=tree][aria-label="Scene"]'],
  ['state tree', '[role=tree][aria-label="State tree"]'],
  ['inspector', '[aria-label="Inspector"], [aria-label="Cell"]'],
  ['sequence band', '[role=tabpanel][aria-label="Sequence"], [data-msp-hover-source]'],
  ['trajectory bar', 'section[aria-label="Playback"]'],
  ['command menu', '[role=dialog]'],
];

/** Where a script ran from: the UI's packages (Prism, Seq*) or Mol*'s (the rest of the workspace and its deps). */
const isUiSource = (url: string) => /\/packages\/(prism|sequences)\//.test(url) && !/\/_stories\//.test(url);

export interface AnimationResult {
  setup: AnimationSetup;
  frames: number;
  frameMs: number;
  p95FrameMs: number;
  longFrames: number;
  /** Long animation frames' script ms per frame, split by where the script came from. */
  uiScriptMsPerFrame: number;
  molstarScriptMsPerFrame: number;
  writesPerFrame: number;
  rendersPerFrame: number;
  mutationsPerFrame: Record<string, number>;
}

const counts = { write: 0, render: 0 };
const frame = () => new Promise<void>((r) => requestAnimationFrame(() => r()));
const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));

/** Puts the viewer's panels in a setup: all closed, then one tab shown or the command menu opened. */
async function arrange(viewer: PrismViewer, shadow: ShadowRoot, setup: AnimationSetup) {
  const shell = viewer.views.shell!;
  shadow
    .querySelector<HTMLElement>('[role=dialog]')
    ?.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, composed: true }));
  shell.closeOverlay();
  for (const d of Docks) viewer.managers.navigation.setDock(d, { collapsed: true });
  await wait(300);
  if (setup === 'structure') shell.show('structure');
  else if (setup === 'state-tree') shell.show('state-tree.tab.tree');
  else if (setup === 'sequence') shell.show('sequence');
  else if (setup === 'menu')
    shadow.querySelector<HTMLElement>('button[title="Command menu"], [aria-label="Command menu"]')?.click();
  await wait(600);
}

/** Measures `seconds` of the playing animation. */
async function measure(shadow: ShadowRoot, setup: AnimationSetup, seconds: number): Promise<AnimationResult> {
  const intervals: number[] = [];
  let ui = 0;
  let molstar = 0;
  let longFrames = 0;
  const mutations: Record<string, number> = {};
  for (const [name] of Areas) mutations[name] = 0;
  mutations.other = 0;
  const observer = new MutationObserver((list) => {
    for (const r of list) {
      const el = r.target.nodeType === Node.ELEMENT_NODE ? (r.target as Element) : r.target.parentElement;
      const area = el ? Areas.find(([, sel]) => el.closest(sel)) : undefined;
      mutations[area ? area[0] : 'other']++;
    }
  });
  observer.observe(shadow, { subtree: true, childList: true, attributes: true, characterData: true });
  const loaf = PerformanceObserver.supportedEntryTypes?.includes('long-animation-frame')
    ? new PerformanceObserver((list) => {
        for (const e of list.getEntries() as any[]) {
          longFrames++;
          for (const s of e.scripts ?? []) {
            if (isUiSource(s.sourceURL ?? '')) ui += s.duration;
            else molstar += s.duration;
          }
        }
      })
    : undefined;
  loaf?.observe({ type: 'long-animation-frame' });
  counts.write = counts.render = 0;
  const end = performance.now() + seconds * 1000;
  let last = performance.now();
  await frame();
  last = performance.now();
  while (performance.now() < end) {
    await frame();
    const now = performance.now();
    intervals.push(now - last);
    last = now;
  }
  const writes = counts.write;
  const renders = counts.render;
  observer.disconnect();
  loaf?.disconnect();
  const n = Math.max(1, intervals.length);
  const sorted = intervals.slice().sort((a, b) => a - b);
  const r = (x: number, d = 2) => +x.toFixed(d);
  const perFrame: Record<string, number> = {};
  for (const [k, v] of Object.entries(mutations)) if (v) perFrame[k] = r(v / n, 1);
  return {
    setup,
    frames: intervals.length,
    frameMs: r(intervals.reduce((a, b) => a + b, 0) / n),
    p95FrameMs: r(sorted[Math.floor(n * 0.95)] ?? 0),
    longFrames,
    uiScriptMsPerFrame: r(ui / n, 3),
    molstarScriptMsPerFrame: r(molstar / n, 3),
    writesPerFrame: r(writes / n, 1),
    rendersPerFrame: r(renders / n, 1),
    mutationsPerFrame: perFrame,
  };
}

/** The real viewer with 1GRM, its models playing at `ModelsPerSecond`, five setups to measure; a readout beside it. */
export const UiCostWhileAnimating: StoryObj = {
  name: 'UI cost while animating',
  parameters: { a11y: { disable: true, test: 'off' }, layout: 'fullscreen' },
  render: () =>
    story(() => {
      // Counted from the first binding the viewer makes (the probe is read when a binding is made).
      RuntimeProbe.on = (kind) => {
        counts[kind]++;
      };
      const container = document.createElement('div');
      container.className = 'msp-root';
      container.style.cssText = 'position: fixed; inset: 0; display: grid; grid-template-columns: 1fr 380px;';
      const host = document.createElement('div');
      host.style.cssText = 'position: relative; min-width: 0; height: 100vh;';
      const panel = document.createElement('div');
      panel.style.cssText =
        'overflow: auto; padding: var(--msp-space-4); font-size: var(--msp-font-size-sm); color: var(--msp-color-text); background: var(--msp-color-surface);';
      container.append(host, panel);
      document.body.append(container);
      let viewer: PrismViewer | undefined;
      let gone = false;
      const status = observableState('Loading 1GRM…');
      const results = observableState<readonly AnimationResult[]>([]);
      const running = observableState(false);
      const ready = (async () => {
        await frame();
        const made = await createPrismViewer({ target: host, spec: DefaultPrismViewerSpec() });
        if (gone) {
          await made.dispose();
          return undefined;
        }
        viewer = made;
        const r = await made.api.load.url!({ url: ExampleFiles.grm, format: 'mmcif' });
        if (!r.ok) status.set(r.message);
        const t = await made.api.trajectory.get!({});
        const count = (t.ok && (t.value as any)?.structures?.[0]?.count) || 1;
        const played = await made.api.trajectory.play!({ mode: 'loop', seconds: Math.max(1, count / ModelsPerSecond) });
        status.set(
          played.ok
            ? `1GRM, ${count} models, playing at ${ModelsPerSecond} models/s.`
            : `1GRM, ${count} models; play: ${played.message}`,
        );
        return made;
      })();
      const run = async (only?: readonly AnimationSetup[], seconds = 5) => {
        const v = await ready;
        if (!v || !host.shadowRoot) return [];
        running.set(true);
        const out: AnimationResult[] = [];
        for (const setup of only ?? Setups) {
          status.set(`Measuring: ${setup}…`);
          await arrange(v, host.shadowRoot, setup);
          out.push(await measure(host.shadowRoot, setup, seconds));
          results.set(out.slice());
        }
        status.set('Done.');
        running.set(false);
        return out;
      };
      (window as any).__mspAnim = {
        ready: () => ready.then((v) => !!v),
        viewer: () => viewer,
        run,
        arrange: async (setup: AnimationSetup) => {
          const v = await ready;
          if (v && host.shadowRoot) await arrange(v, host.shadowRoot, setup);
        },
        measure: async (setup: AnimationSetup, seconds = 5) =>
          host.shadowRoot && measure(host.shadowRoot, setup, seconds),
      };
      const disposeReadout = render(
        () => (
          <div style="display: grid; gap: var(--msp-space-3);">
            <strong>UI cost while animating</strong>
            <span>{status}</span>
            <button
              type="button"
              style="justify-self: start; height: 28px; padding: 0 var(--msp-space-4); border-radius: 6px; background: var(--msp-color-control);"
              disabled={running}
              onClick={() => void run()}
            >
              {derive(running.asReadonly, (r) => (r ? 'Running…' : 'Run (5 s a setup)'))}
            </button>
            <table style="border-collapse: collapse; font-variant-numeric: tabular-nums;">
              <thead>
                <tr>
                  <th style="text-align: left">setup</th>
                  <th>frame ms</th>
                  <th>UI script</th>
                  <th>Mol* script</th>
                  <th>writes</th>
                  <th>renders</th>
                </tr>
              </thead>
              <tbody>
                <For each={results}>
                  {(r) => (
                    <tr title={JSON.stringify(r.mutationsPerFrame)}>
                      <td>{r.setup}</td>
                      <td style="text-align: right">{`${r.frameMs} (p95 ${r.p95FrameMs})`}</td>
                      <td style="text-align: right">{r.uiScriptMsPerFrame}</td>
                      <td style="text-align: right">{r.molstarScriptMsPerFrame}</td>
                      <td style="text-align: right">{r.writesPerFrame}</td>
                      <td style="text-align: right">{r.rendersPerFrame}</td>
                    </tr>
                  )}
                </For>
              </tbody>
            </table>
            <span style="color: var(--msp-color-text-subtle)">
              Per frame. Script: long animation frames' attribution (frames over 50 ms only); hover a row for its DOM
              mutations by area.
            </span>
          </div>
        ),
        panel,
      );
      onCleanup(() => {
        gone = true;
        RuntimeProbe.on = undefined;
        disposeReadout();
        void viewer?.dispose();
        container.remove();
        (window as any).__mspAnim = undefined;
      });
      return (
        <p style="padding: var(--msp-space-4)">The viewer and its readout are mounted outside Storybook's root.</p>
      );
    }),
};
