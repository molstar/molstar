/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { For, Legend } from '@molstar/prism-ui';
import type { SequenceViewer } from './model.js';

import { derive } from '@molstar/prism-ui/reactive/state';
/** The legend of the sequence's colors, in Prism's `Legend`: the card the Sequence menu's "Legend…" opens. */
export function LegendCard(props: { viewer: SequenceViewer }) {
  const { legend } = props.viewer;
  return (
    <div style="display: grid; gap: var(--msp-space-4); padding: var(--msp-space-3)">
      <Legend
        title={derive(legend, (l) => l.title)}
        items={derive(legend, (l) => l.items)}
        ramps={derive(legend, (l) => l.ramps)}
      />
      <For each={derive(legend, (l) => l.sections ?? [])}>
        {(section) => <Legend title={section.title} items={section.items} />}
      </For>
    </div>
  );
}
