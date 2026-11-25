/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ExampleView } from '@molstar/seqx-examples';
import { mountLegend } from '@molstar/seqx-viewer/mount-legend';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { followSvg } from '../../follow-svg.js';
import { composeViewer } from '../selection/compose.js';
import { Channels, cycleStyle, HighlightConfig, describeChannels, toggleChannel } from './compose.js';

// The highlights story's root: p53 with three highlight channels that overlap (a domain, a motif search result, a host's
// own), each set and cleared from a button and given its next style (tint, band, outline) from another. Where channels
// overlap the band splits into equal stripes in the order they were set; one alone takes the whole band. The Legend
// part lists the active channels with their colors, each removable; the SVG beside the canvas draws the same channels.

const Box = 'display: grid; gap: var(--msui-space-2); min-width: 0; align-content: start;';
const Title =
  'margin: 0; font-size: var(--msui-text-sm); font-weight: var(--msui-weight-strong); color: var(--msui-text-muted);';

function section(title: string, body: HTMLElement): HTMLElement {
  const el = document.createElement('div');
  el.style.cssText = Box;
  const h = document.createElement('h4');
  h.style.cssText = Title;
  h.textContent = title;
  el.append(h, body);
  return el;
}

class HighlightsRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-highlights-root';

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; align-content: start; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent =
      'Set the three channels one by one: where two or more overlap, the band is split into equal stripes in the order they were set; where one lights alone, it takes the whole band. Each channel has a style: a tint, a band along the foot, or an outline. The legend lists the active channels, each with a clear.';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const svgHost = document.createElement('div');
    svgHost.style.cssText = 'min-width: 0; overflow-x: auto;';
    const legendHost = document.createElement('div');
    const status = document.createElement('p');
    status.style.cssText =
      'margin: 0; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); color: var(--msui-text-muted);';
    const tools = document.createElement('div');
    tools.style.cssText = 'display: flex; flex-wrap: wrap; gap: var(--msui-space-2); align-items: baseline;';
    const pair = document.createElement('div');
    pair.style.cssText =
      'display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: var(--msui-space-5);';
    pair.append(section('Canvas', canvasHost), section('SVG', svgHost));
    root.append(hint, tools, pair, section('Legend', legendHost), status);

    const viewer = composeViewer(view, { target: canvasHost, noRails: true, config: HighlightConfig });
    onCleanup(() => viewer.retire());
    for (const c of Channels) {
      tools.append(
        Button({ children: `Toggle ${c.title}`, onClick: () => void toggleChannel(viewer, c) }) as HTMLElement,
        Button({ children: `${c.name}: next style`, onClick: () => void cycleStyle(viewer, c) }) as HTMLElement,
      );
    }
    onCleanup(mountLegend(legendHost, viewer));

    const svg = followSvg(viewer, svgHost, { prefix: 'hl', also: [viewer.highlight.version] });
    onCleanup(() => svg.dispose());
    const showStatus = () => {
      status.textContent = describeChannels(viewer);
    };
    onCleanup(viewer.highlight.version.listen(showStatus));
    void viewer
      .load(view.build())
      .then(() => viewer.settle())
      .then(() => {
        showStatus();
        return svg.draw();
      })
      .catch((e: unknown) => {
        status.textContent = `The view did not load: ${e instanceof Error ? e.message : String(e)}`;
      });
    return root;
  }
}

defineElements([HighlightsRoot]);

/** A Storybook `render` result: the view with the buttons that set its highlight channels. */
export function highlightsStory(view: ExampleView): HTMLElement {
  const el = document.createElement(HighlightsRoot.tag) as HighlightsRoot;
  el.model = view;
  return el;
}
