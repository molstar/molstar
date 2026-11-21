/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { modifiersOf } from '@molstar/seqx-viewer/bindings';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { requestFrameCallback } from '@molstar/ui/dom/frame';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { HostSection, hostDrawnSpec, hostPicture } from './compose.js';

// The host-drawn navigator story (the viewer's _spec "Models", NavigatorModel, the contract): a navigator the host
// draws itself in plain DOM, no canvas: the viewer's `navigator(def)` model holds the rules (the window's geometry,
// floor, gestures, keys, idle rule, readout, sliders' values); the host maps its pointer and keys to the model's methods
// and draws the domains and the windows as elements from the model's last frame. The viewer ends with the story.

export interface HostNavigatorStory {
  readonly hint: string;
}

/** The navigator's height (CSS px): the host's choice. */
const Height = 36;

class HostNavigatorRoot extends UiElement<HostNavigatorStory> {
  static readonly tag = 'seqx-host-navigator-root';

  protected render(story: HostNavigatorStory): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; align-content: start; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent = story.hint;
    // the host's own element: one tab stop, a slider whose value the model says
    const nav = document.createElement('div');
    nav.tabIndex = 0;
    nav.setAttribute('role', 'slider');
    nav.setAttribute('aria-label', 'Host-drawn navigator');
    nav.style.cssText = `position: relative; height: ${Height}px; overflow: hidden; touch-action: none; border-radius: var(--msui-radius-sm, 4px); background: var(--msui-surface-sunken, var(--msui-surface)); outline-offset: 2px;`;
    const readout = document.createElement('div');
    readout.style.cssText = 'min-height: 1.5em; color: var(--msui-text-muted); font-variant-numeric: tabular-nums;';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    root.append(hint, nav, readout, canvasHost);

    const viewer: SeqViewer = createSeqViewer(DefaultSeqViewerSpec(), canvasHost);
    const model = viewer.navigator({ id: 'host', of: HostSection, content: 'axis' });
    let frame: number | undefined;
    const paint = () => {
      frame = undefined;
      model.resize(nav.clientWidth, Height);
      if (!model.draw()) return;
      const pic = hostPicture(model);
      const els: HTMLElement[] = [];
      for (const b of pic.blocks) {
        const el = document.createElement('div');
        el.textContent = b.label;
        el.style.cssText = `position: absolute; left: ${b.x}px; width: ${b.width}px; top: 8px; bottom: 8px; border-radius: 3px; background: var(--msui-accent-soft, color-mix(in srgb, var(--msui-accent) 25%, transparent)); color: var(--msui-text); font-size: 10px; line-height: ${Height - 16}px; text-align: center; overflow: hidden; pointer-events: none;`;
        els.push(el);
      }
      for (const w of pic.windows) {
        const el = document.createElement('div');
        el.style.cssText = `position: absolute; left: ${w.x}px; width: ${w.width}px; top: 1px; bottom: 1px; box-sizing: border-box; border: 2px solid ${w.idle ? 'var(--msui-line)' : 'var(--msui-accent)'}; border-radius: 4px; background: color-mix(in srgb, var(--msui-accent) 10%, transparent); pointer-events: none;`;
        els.push(el);
      }
      nav.replaceChildren(...els);
      nav.style.cursor = model.cursor;
      readout.textContent = model.readout.value ?? 'Hover the navigator';
    };
    const schedule = () => {
      if (frame === undefined) frame = requestFrameCallback(paint);
    };
    const say = () => nav.setAttribute('aria-valuetext', model.values.value[model.focused.value]?.text ?? '');
    const offs = [model.version.listen(schedule), viewer.drawn.listen(schedule), model.values.listen(say)];

    // the host's events, each one of the model's methods
    const at = (e: MouseEvent) => {
      const r = nav.getBoundingClientRect();
      return [e.clientX - r.left, e.clientY - r.top] as const;
    };
    const touch = (e: PointerEvent) => ({ touch: e.pointerType === 'touch' });
    nav.addEventListener('pointermove', (e) => model.move(...at(e), touch(e)));
    nav.addEventListener('pointerleave', () => model.leave());
    nav.addEventListener('pointerdown', (e) => {
      if (e.button !== 0) return;
      if (model.down(...at(e), touch(e))) nav.setPointerCapture(e.pointerId);
      nav.focus({ preventScroll: true });
    });
    nav.addEventListener('pointerup', () => model.up());
    nav.addEventListener('pointercancel', () => model.cancel());
    nav.addEventListener('click', (e) => model.click(...at(e), e.detail));
    nav.addEventListener('dblclick', (e) => model.doubleClick(...at(e)));
    nav.addEventListener('focus', () => model.setFocus(true));
    nav.addEventListener('blur', () => model.setFocus(false));
    nav.addEventListener(
      'wheel',
      (e) => {
        if (model.wheel(...at(e), e.deltaX, e.deltaY, modifiersOf(e, viewer.input.mac))) e.preventDefault();
      },
      { passive: false },
    );
    nav.addEventListener('keydown', (e) => {
      if (model.key(e.key, e)) e.preventDefault();
    });
    const resized = new ResizeObserver(schedule);
    resized.observe(nav);
    onCleanup(() => {
      for (const off of offs) off();
      resized.disconnect();
      viewer.retire();
    });
    void viewer.load(hostDrawnSpec()).then(schedule);
    return root;
  }
}

defineElements([HostNavigatorRoot]);

/** A Storybook `render` result: the host's navigator above the viewer it moves. */
export function hostNavigatorStory(story: HostNavigatorStory): HTMLElement {
  const el = document.createElement(HostNavigatorRoot.tag) as HostNavigatorRoot;
  el.model = story;
  return el;
}
