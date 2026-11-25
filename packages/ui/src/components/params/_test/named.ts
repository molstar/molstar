/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** The field a row's label names (`aria-labelledby="<labelId> <labelId>-value"`, or an `aria-label`). */
export function fieldNamed(container: ParentNode, label: string): HTMLButtonElement | null {
  for (const b of container.querySelectorAll<HTMLButtonElement>('button[aria-labelledby], button[aria-label]')) {
    const by = b.getAttribute('aria-labelledby');
    const name = by ? document.getElementById(by.split(/\s+/)[0])?.textContent?.trim() : b.getAttribute('aria-label');
    if (name === label) return b;
  }
  return null;
}
