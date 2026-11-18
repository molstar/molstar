// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { KeyPlatform } from '@molstar/ui/commands/keys';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { RecordingContext } from '../canvas/_test/recording.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { DefaultToolbar, toolbarButtons } from '../toolbar.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { ProteinSection as Protein, proteinOnGene } from './fixture.js';
import { type Frames, installFrames } from './frames.js';

// The default toolbar (_spec/README.md "Input", Toolbar): the buttons as the DOM-free model gives them, and the bar they
// make in a DOM: Back and Forward with their keys in the tooltip, disabled while their side of the view's history is
// empty, running the commands.

let frames: Frames;
let viewer: SeqViewer | undefined;
let mac: boolean;

beforeEach(() => {
  frames = installFrames();
  mac = KeyPlatform.mac;
  KeyPlatform.mac = false;
  vi.stubGlobal(
    'ResizeObserver',
    class {
      observe() {}
      unobserve() {}
      disconnect() {}
    },
  );
  Object.defineProperty(document.defaultView!, 'matchMedia', {
    configurable: true,
    writable: true,
    value: () => ({ matches: false, addEventListener() {}, removeEventListener() {} }),
  });
  vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as any).__ctx ??= new RecordingContext());
  });
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  KeyPlatform.mac = mac;
  frames.restore();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
  Reflect.deleteProperty(document.defaultView!, 'matchMedia');
  document.body.innerHTML = '';
});

async function make(mounted = false) {
  const target = document.createElement('div');
  document.body.appendChild(target);
  viewer = createSeqViewer(
    { ...DefaultSeqViewerSpec(), config: [[viewerItem('viewport.maxPxPerPosition'), 1000]] },
    mounted ? target : undefined,
  );
  await viewer.load(proteinOnGene());
  viewer.resize(600);
  await viewer.settle();
  frames.flush();
  return { target, viewer };
}

describe('the toolbar’s buttons, headless', () => {
  it('are disabled with a reason while their side of the history is empty, across back and forward', async () => {
    const { viewer: v } = await make();
    const [back, forward] = toolbarButtons(v);
    expect(DefaultToolbar).toEqual(['view.back', 'view.forward']);
    expect([back.command, forward.command]).toEqual(['view.back', 'view.forward']);
    expect([back.disabled.value, forward.disabled.value]).toEqual(['No view to go back to', 'No view to go forward to']);
    const told: unknown[] = [];
    back.disabled.listen((d) => told.push(d));
    v.execute('view.focus', { section: Protein, focus: { start: 1, end: 3 } });
    expect([back.disabled.value, forward.disabled.value]).toEqual([false, 'No view to go forward to']);
    expect(await back.run()).toBe('Back');
    expect([back.disabled.value, forward.disabled.value]).toEqual(['No view to go back to', false]);
    expect(await forward.run()).toBe('Forward');
    expect([back.disabled.value, forward.disabled.value]).toEqual([false, 'No view to go forward to']);
    // told once a change of its own, never for what it did not change
    expect(told).toEqual([false, 'No view to go back to', false]);
  });

  it('say what they do and their key: the first of `view.keys`, as the platform writes it', async () => {
    const { viewer: v } = await make();
    const [back, forward] = toolbarButtons(v);
    expect([back.title, forward.title]).toEqual(['Back (Alt+←)', 'Forward (Alt+→)']);
    KeyPlatform.mac = true;
    expect(toolbarButtons(v)[0].title).toBe('Back (⌥←)');
    KeyPlatform.mac = false;
    // a command the bar does not know is left out
    expect(toolbarButtons(v, ['view.fit', 'view.back']).map((b) => b.command)).toEqual(['view.back']);
  });

  it('show no key for a command whose keys the application removed', async () => {
    const v = createSeqViewer({
      ...DefaultSeqViewerSpec(),
      ui: undefined,
      config: [[viewerItem('view.keys'), { 'view.back': [] }]],
    });
    viewer = v;
    expect(toolbarButtons(v)[0]).toMatchObject({ title: 'Back', key: undefined });
  });
});

describe('the mounted toolbar', () => {
  const buttons = (target: HTMLElement) => [...target.querySelectorAll<HTMLButtonElement>('[data-toolbar-button]')];

  it('is above the canvas, its buttons disabled with the reason until there is a step, and run the commands', async () => {
    const { target, viewer: v } = await make(true);
    expect(target.firstElementChild?.querySelector('[role="toolbar"]')).not.toBeNull();
    expect(target.querySelector('[role="application"]')).not.toBeNull();
    const [back, forward] = buttons(target);
    // the accessible name is the title, key included
    expect([back.getAttribute('aria-label'), forward.getAttribute('aria-label')]).toEqual([
      'Back (Alt+←)',
      'Forward (Alt+→)',
    ]);
    expect(back.getAttribute('aria-disabled')).toBe('true');
    expect(forward.getAttribute('aria-disabled')).toBe('true');
    expect(back.getAttribute('data-tip') ?? back.title).toContain('No view to go back to');
    v.execute('view.focus', { section: Protein, focus: { start: 1, end: 3 } });
    frames.flush();
    expect(back.getAttribute('aria-disabled')).not.toBe('true');
    expect(back.getAttribute('data-tip') ?? back.title).toContain('Alt+←');
    const focused = [v.view.get(Protein)!.start, v.view.get(Protein)!.end];
    back.click();
    frames.flush();
    expect([v.view.get(Protein)!.start, v.view.get(Protein)!.end]).not.toEqual(focused);
    expect(back.getAttribute('aria-disabled')).toBe('true');
    expect(forward.getAttribute('aria-disabled')).not.toBe('true');
    forward.click();
    expect([v.view.get(Protein)!.start, v.view.get(Protein)!.end]).toEqual(focused);
  });

  it('is left out with no `ui.toolbar`, and retire takes it away with the canvas', async () => {
    const target = document.createElement('div');
    document.body.appendChild(target);
    const bare = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: { canvas: true } }, target);
    expect(target.querySelector('[role="toolbar"]')).toBeNull();
    bare.retire();
    const { target: t2, viewer: v } = await make(true);
    expect(t2.querySelector('[role="toolbar"]')).not.toBeNull();
    v.retire();
    expect(t2.childElementCount).toBe(0);
  });
});
