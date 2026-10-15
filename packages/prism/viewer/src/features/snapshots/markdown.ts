/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ColorLists } from '@molstar/core/util/color/lists';
import { getColorGradient, getColorGradientBanded, parseColorList } from '@molstar/core/util/color/utils';
import type { MarkdownHost } from '@molstar/prism-ui';
import type { ViewerApi } from '../../api/area.js';
import type { MarkdownSource } from './source.js';

// The descriptions' markdown host (B1): Mol*'s markdown extensions behind prism-ui's `Markdown`, as v5's
// `controls/markdown.tsx` puts them behind react-markdown. Command links run through `snapshot.link` (the person's
// API); images come from the snapshot area's markdown source (the audio player has its own bar in the card). The
// inline renderers are v5's defaults (color swatch, color palette); an extension's own `reactRenderFn` is React and is
// not shown.

export function createMarkdownHost(api: ViewerApi, source: MarkdownSource): MarkdownHost {
  return {
    link(href) {
      if (!source.command(href)) return undefined;
      return (event) => void api.snapshot?.link({ href, event });
    },
    image(src, alt) {
      const args = source.command(src);
      if (args) return renderCommand(args) ?? document.createTextNode(alt);
      return source.resolve(src);
    },
    // The audio player is the card's own bar (views/state/audio-bar.tsx), not the browser's player in the text.
  };
}

/** v5's `DefaultMarkdownExtensionRenderers`: `![](!color-swatch=red)`, `![](!color-palette-name=viridis)`. */
function renderCommand(args: Readonly<Record<string, string>>): Node | undefined {
  const swatch = args['color-swatch'];
  if (swatch) {
    return span({
      display: 'inline-block',
      width: '0.75em',
      height: '0.75em',
      borderRadius: '25%',
      backgroundColor: swatch,
    });
  }
  const name = args['color-palette-name'];
  const colors = args['color-palette-colors'];
  if (!name && !colors) return undefined;
  const list = colors ? parseColorList(colors) : ColorLists[name.toLowerCase() as keyof typeof ColorLists]?.list;
  if (!list?.length) return undefined;
  return span({
    display: 'inline-block',
    minWidth: args['color-palette-width'] ?? '150px',
    height: args['color-palette-height'] ?? '0.5em',
    borderRadius: '2px',
    background: ('color-palette-discrete' in args ? getColorGradientBanded : getColorGradient)(list),
  });
}

function span(style: Partial<CSSStyleDeclaration>): HTMLSpanElement {
  const el = document.createElement('span');
  Object.assign(el.style, style);
  return el;
}
