/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '@molstar/ui/components/button/button';
import { lucideEye } from '@molstar/ui/components/icon/_generated/lucide/eye';
import { lucideEyeOff } from '@molstar/ui/components/icon/_generated/lucide/eye-off';
import { installTooltips } from '@molstar/ui/components/tooltip/tooltip';
import { render } from '@molstar/ui/dom/render';
import { onCleanup } from '@molstar/ui/dom/scope';

// The Legend part's eyes (_spec/README.md "Legends", A layer's eye): `@molstar/ui`'s ghost icon button and its tooltip,
// a chunk of their own that the part loads the first time it draws an eye (the default viewer's first chunk carries
// none of the button, the tooltip or the DOM runtime they need). Each returns what takes it away.

/** A Lucide icon's markup as an element (one icon imported by itself, never the registry of every icon). */
function iconOf(doc: Document, markup: string): SVGSVGElement {
  const svg = doc.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', '0 0 24 24');
  svg.setAttribute('fill', 'none');
  svg.setAttribute('stroke', 'currentColor');
  svg.setAttribute('stroke-width', '2');
  svg.setAttribute('stroke-linecap', 'round');
  svg.setAttribute('stroke-linejoin', 'round');
  svg.setAttribute('aria-hidden', 'true');
  svg.innerHTML = markup;
  return svg;
}

export interface EyeButton {
  /** Its name: its tooltip and accessible name ("Hide Contacts"). */
  readonly title: string;
  /** The eye-off (its layer hidden), else the eye. */
  readonly off: boolean;
  /** It opens a menu (a legend over several layers). */
  readonly menu: boolean;
  onClick(button: HTMLElement): void;
}

/** Renders an eye into `holder`. */
export function renderEye(holder: HTMLElement, eye: EyeButton): () => void {
  return render(() => {
    const b = Button({
      variant: 'ghost',
      size: 'xs',
      icon: iconOf(holder.ownerDocument, eye.off ? lucideEyeOff : lucideEye),
      title: eye.title,
      onClick: (ev) => eye.onClick(ev.currentTarget as HTMLElement),
    }) as HTMLElement;
    if (eye.menu) b.setAttribute('aria-haspopup', 'menu');
    return b;
  }, holder);
}

/** The tooltip the eyes' `data-tip` shows, in `root`. */
export function eyeTooltips(root: HTMLElement): () => void {
  return render(() => {
    onCleanup(installTooltips(root));
    return null;
  }, root);
}
