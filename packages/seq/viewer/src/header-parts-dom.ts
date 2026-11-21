/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '@molstar/ui/components/button/button';
import { Icon } from '@molstar/ui/components/icon/icon';
import { render } from '@molstar/ui/dom/render';
import type {
  HeaderMenuItem,
  HeaderPartsHeader,
  HeaderPartsPortFactory,
  HeaderPartsState,
  HeaderTarget,
} from './header-entries.js';

// The headers' reference DOM parts (_spec/README.md "Ports", Header parts): `@molstar/ui` buttons over the boxes the
// layout gives (`viewer.headers`), in the canvas's element: the action slots and the "⋯" of a header whose actions show
// (a hover, the cursor; so a 200-track view holds a few buttons, not 600), and a button over each chip the canvas draws.
// A press runs the slot's entry as the canvas's click does (`HeadersBehavior.run`); the "⋯" and a slot with a submenu
// open the menu at it. Keyboard: the cursor on a header moves the focus into its parts (one tab stop, roving: ← → between
// them, Enter or Space runs), Esc gives it back to the canvas, the cursor where it was. Tab is never captured.

const at = (el: HTMLElement, b: { x: number; y: number; width: number; height: number }) =>
  Object.assign(el.style, {
    position: 'absolute',
    left: `${b.x}px`,
    top: `${b.y}px`,
    width: `${b.width}px`,
    height: `${b.height}px`,
    minWidth: '0',
    padding: '0',
  });

const sameTarget = (a: HeaderTarget, b: HeaderTarget) => a.section === b.section && a.track === b.track;

/** @internal `@molstar/ui` buttons over the headers' boxes (the lazy chunk `headerPartsDom()` loads). */
export function partsDomPort(): HeaderPartsPortFactory {
  return (container, behavior) => {
    const doc = container.ownerDocument;
    const layer = doc.createElement('div');
    layer.className = 'seqx-header-parts';
    Object.assign(layer.style, { position: 'absolute', left: '0', top: '0', width: '0', height: '0' });
    container.append(layer);
    // a part's own events are its own: the canvas below sees no press, click, move or key there (stopped at the parts'
    // root, where `@molstar/ui` delegates the buttons' handlers)
    const own = (e: Event) => e.stopPropagation();
    const types = ['pointerdown', 'pointerup', 'click', 'dblclick', 'pointermove', 'contextmenu', 'keydown', 'wheel'];
    for (const type of types) layer.addEventListener(type, own);
    let dispose: (() => void) | undefined;
    let last: HeaderPartsState = { headers: [] };
    const press = (
      h: HeaderPartsHeader,
      it: HeaderMenuItem | undefined,
      box: { x: number; y: number; height: number },
    ) => {
      const place = { x: box.x, y: box.y + box.height };
      if (!it) return void behavior.openMenu(h.geometry.target, place);
      if (it.items && it.kind === 'action') return void behavior.openMenu(h.geometry.target, place, it.name);
      void behavior.run(it.name, h.geometry.target);
    };
    const draw = (state: HeaderPartsState) => {
      dispose?.();
      layer.textContent = '';
      dispose = render(() => {
        const out = doc.createElement('div');
        for (const h of state.headers) {
          const group = doc.createElement('div');
          group.setAttribute('role', 'group');
          group.setAttribute('aria-label', h.geometry.target.track ?? h.geometry.target.section);
          group.dataset.section = h.geometry.target.section;
          if (h.geometry.target.track) group.dataset.track = h.geometry.target.track;
          const buttons: HTMLElement[] = [];
          // the slots' backdrop in the page's color over the label's tail, with the ellipsis the canvas draws when it cuts
          // the label (the text's end under the slots): no overprinting
          if (h.actions && h.geometry.slots.length) {
            const box = h.geometry.slots[0].box;
            const x0 = Math.min(...h.geometry.slots.map((s) => s.box.x)) - 2;
            const end = Math.max(...h.geometry.slots.map((s) => s.box.x + s.box.width));
            const cut = h.geometry.labelEnd !== undefined && h.geometry.labelEnd > x0;
            const grow = cut ? (h.geometry.ellipsis ?? 0) : 0;
            const back = doc.createElement('div');
            back.setAttribute('aria-hidden', 'true');
            back.dataset.backdrop = '';
            at(back, { x: x0 - grow, y: box.y, width: end - x0 + grow, height: box.height });
            Object.assign(back.style, {
              background: 'var(--msui-bg)',
              pointerEvents: 'none',
              color: 'var(--msui-text-muted)',
              font: 'inherit',
              fontSize: '12px',
              lineHeight: `${box.height}px`,
              overflow: 'hidden',
            });
            if (cut) back.textContent = '…';
            group.append(back);
          }
          if (h.actions)
            for (const s of h.geometry.slots) {
              const it = h.items.find((i) => i.name === s.name);
              const b = Button({
                variant: 'ghost',
                size: 'sm',
                icon: Icon({ name: it?.icon ?? 'ellipsis' }),
                title: it ? it.label : 'Actions',
                ...(it && !it.enabled ? { 'aria-description': it.why } : {}),
                onClick: () => press(h, it, s.box),
              }) as HTMLElement;
              b.dataset.part = s.name;
              if (it && !it.enabled) b.setAttribute('aria-disabled', 'true');
              at(b, s.box);
              buttons.push(b);
            }
          for (const s of h.geometry.states) {
            const it = h.items.find((i) => i.name === s.name);
            const b = doc.createElement('button');
            b.type = 'button';
            b.dataset.part = s.name;
            b.setAttribute('aria-label', it?.command ? `${s.text}: ${it.label}` : s.text);
            // the chip itself (the canvas draws nothing in its box): its words, or a mark's squircle where the column
            // leaves no room for them (an export draws the canvas's)
            if (!s.mark) b.textContent = s.text;
            Object.assign(b.style, {
              background: s.mark ? 'var(--msui-accent)' : 'transparent',
              color: 'var(--msui-accent)',
              border: '0',
              borderRadius: s.mark ? '2px' : '0',
              padding: '0',
              font: 'inherit',
              fontSize: '11px',
              whiteSpace: 'nowrap',
              cursor: it?.command ? 'pointer' : 'default',
            });
            at(b, s.box);
            b.addEventListener('click', () => it?.command && press(h, it, s.box));
            buttons.push(b);
          }
          // one tab stop a header, roving
          buttons.forEach((b, i) => {
            b.tabIndex = i === 0 ? 0 : -1;
            group.append(b);
          });
          group.addEventListener('keydown', (e) => {
            const i = buttons.indexOf(e.target as HTMLElement);
            if (e.key === 'ArrowRight' || e.key === 'ArrowLeft') {
              const next = buttons[(i + (e.key === 'ArrowRight' ? 1 : buttons.length - 1)) % buttons.length];
              for (const b of buttons) b.tabIndex = b === next ? 0 : -1;
              next?.focus();
              e.preventDefault();
            } else if (e.key === 'Escape') {
              container.focus({ preventScroll: true });
              e.preventDefault();
              e.stopPropagation();
            }
          });
          out.append(group);
        }
        return out;
      }, layer);
    };
    return {
      update(state) {
        last = state;
        draw(state);
      },
      focus(target) {
        const h = last.headers.find((x) => sameTarget(x.geometry.target, target));
        if (!h) return;
        const group = [...layer.querySelectorAll<HTMLElement>('[role="group"]')].find(
          (g) => g.dataset.section === target.section && g.dataset.track === target.track,
        );
        group?.querySelector<HTMLElement>('button')?.focus({ preventScroll: true });
      },
      dispose() {
        dispose?.();
        for (const type of types) layer.removeEventListener(type, own);
        layer.remove();
      },
    };
  };
}
