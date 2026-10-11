/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineElements, MspElement } from '../dom/element.js';
import type { Child } from '../dom/jsx-types.js';

/** Hosts a story's JSX in an element, so its bindings are disposed when Storybook removes it. */
class MspStoryRoot extends MspElement<() => Child> {
  static tag = 'msp-story-root';
  protected render(content: () => Child) {
    return content();
  }
}

defineElements([MspStoryRoot]);

/** A Storybook `render` result for JSX content. */
export function story(content: () => Child): HTMLElement {
  const el = document.createElement(MspStoryRoot.tag) as MspStoryRoot;
  el.model = content;
  return el;
}
