/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ExampleView, P53Views } from '@molstar/seqx-examples';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { followSvg } from '../../follow-svg.js';
import {
  cameraText,
  composeStory,
  endScrub,
  exportSteps,
  nextStep,
  play,
  scrub,
  scrubTo,
  storyRows,
} from './story-compose.js';

// p53 as a story: each step a view and a small change of the spec, played by the page's clock (the default adapter of
// the animation port), scrubbed by a slider (the application's clock), and followed by a stand-in structure camera
// driven by the same plans. The SVG beside the canvas follows what is shown; Export draws an SVG per step headlessly.

const Row = 'display: flex; gap: var(--msui-space-2); align-items: center; flex-wrap: wrap;';

class StoryRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-p53-story-root';

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; align-content: start; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent =
      'Play a step, or Next to walk the story; drag the slider to scrub the transition into the step chosen. The camera beside the canvas follows the same plans: a structure viewer would move with it.';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const svgHost = document.createElement('div');
    svgHost.style.cssText = 'min-width: 0; overflow-x: auto;';
    const status = document.createElement('p');
    status.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    const say = (m: string) => {
      status.textContent = m;
    };

    const s = composeStory(view, { target: canvasHost });
    onCleanup(() => s.dispose());

    // the stand-in camera: a dial whose hand turns to where the protein's centre in view is, longer as it zooms in
    const camera = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    camera.setAttribute('viewBox', '-50 -50 100 100');
    camera.setAttribute('width', '96');
    camera.setAttribute('height', '96');
    camera.setAttribute('role', 'img');
    const dial = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
    dial.setAttribute('r', '44');
    dial.setAttribute('fill', 'none');
    dial.setAttribute('stroke', 'var(--msui-line)');
    const hand = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    hand.setAttribute('stroke', 'var(--msui-accent, currentColor)');
    hand.setAttribute('stroke-width', '3');
    hand.setAttribute('stroke-linecap', 'round');
    camera.append(dial, hand);
    const cameraLabel = document.createElement('span');
    cameraLabel.style.cssText = 'color: var(--msui-text-muted);';
    const drawCamera = () => {
      const r = Math.min(44, 12 + 8 * Math.log2(Math.max(1, s.camera.zoom)));
      const a = ((s.camera.angle - 90) * Math.PI) / 180;
      hand.setAttribute('x2', String(r * Math.cos(a)));
      hand.setAttribute('y2', String(r * Math.sin(a)));
      camera.setAttribute('aria-label', cameraText(s.camera));
      cameraLabel.textContent = cameraText(s.camera);
    };
    onCleanup(s.viewer.animation.drawn.listen(drawCamera));

    const steps = document.createElement('div');
    steps.style.cssText = Row;
    const slider = document.createElement('input');
    slider.type = 'range';
    slider.min = '0';
    slider.max = '1000';
    slider.value = '0';
    slider.setAttribute('aria-label', 'The transition into the step chosen, scrubbed');
    const scrubTarget = document.createElement('select');
    scrubTarget.setAttribute('aria-label', 'The step to scrub to');
    const fill = () => {
      steps.replaceChildren();
      scrubTarget.replaceChildren();
      for (const { id, label } of storyRows(s)) {
        steps.append(Button({ size: 'sm', children: label, onClick: () => say(play(s, id)) }) as HTMLElement);
        const o = document.createElement('option');
        o.value = id;
        o.textContent = label;
        scrubTarget.append(o);
      }
      steps.append(
        Button({
          size: 'sm',
          children: 'Next',
          onClick: () => {
            const id = nextStep(s);
            say(id ? play(s, id) : 'The story’s end');
          },
        }) as HTMLElement,
        Button({ size: 'sm', children: 'Leave the story', onClick: () => say(play(s, undefined)) }) as HTMLElement,
      );
    };
    slider.addEventListener('pointerdown', () => {
      if (!s.scrubbing) say(scrubTo(s, scrubTarget.value));
    });
    slider.addEventListener('input', () => {
      if (!s.scrubbing) say(scrubTo(s, scrubTarget.value));
      scrub(s, Number(slider.value) / 1000);
      drawCamera();
    });
    slider.addEventListener('change', () => {
      endScrub(s);
      slider.value = '0';
    });
    const scrubRow = document.createElement('label');
    scrubRow.style.cssText = Row;
    scrubRow.append('Scrub to', scrubTarget, slider);

    const exported = document.createElement('div');
    exported.style.cssText =
      'display: grid; gap: var(--msui-space-3); grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));';
    const exportButton = Button({
      size: 'sm',
      children: 'Export an SVG per step',
      onClick: () => {
        say('Exporting…');
        void exportSteps(view)
          .then((svgs) => {
            exported.replaceChildren();
            for (const { id, svg } of svgs) {
              const box = document.createElement('figure');
              box.style.cssText = 'margin: 0; min-width: 0; overflow-x: auto;';
              const caption = document.createElement('figcaption');
              caption.textContent = id;
              const holder = document.createElement('div');
              // our own output: every string from data or the spec is escaped by the backend
              holder.innerHTML = svg;
              box.append(caption, holder);
              exported.append(box);
            }
            say(`Exported ${svgs.length} steps`);
          })
          .catch((e: unknown) => say(`The export failed: ${e instanceof Error ? e.message : String(e)}`));
      },
    }) as HTMLElement;

    const cameraRow = document.createElement('div');
    cameraRow.style.cssText = Row;
    cameraRow.append(camera, cameraLabel);
    root.append(hint, canvasHost, steps, scrubRow, cameraRow, status, svgHost, exportButton, exported);

    const svg = followSvg(s.viewer, svgHost, { prefix: 'p53s' });
    onCleanup(() => svg.dispose());
    void s.viewer
      .load(view.build())
      .then(() => s.viewer.settle())
      .then(() => {
        fill();
        drawCamera();
        return svg.draw();
      })
      .catch((e: unknown) => say(`The view did not load: ${e instanceof Error ? e.message : String(e)}`));
    return root;
  }
}

defineElements([StoryRoot]);

const meta: Meta = { title: 'Features/Story' };
export default meta;

/** Play the p53 story, scrub a transition with the slider, watch the stand-in camera follow; export a step's SVG. */
export const Story: StoryObj = {
  render: () => {
    const el = document.createElement(StoryRoot.tag) as StoryRoot;
    el.model = P53Views.find((v) => v.id === 'story')!;
    return el;
  },
};
