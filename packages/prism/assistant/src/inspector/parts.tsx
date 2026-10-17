/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { BudgetView, TokensView } from '@molstar/ai/trace';
import { credits, formatUsage, k } from '@molstar/ai/trace';
import { ownsKey } from '@molstar/prism-ui/dom/events';
import { on } from '@molstar/prism-ui/dom/scoped';
import { type Child, type MaybeReactive, Show } from '@molstar/prism-ui';
import { distinctUntilChanged, map, type Observable } from 'rxjs';
import { classes } from './_generated/parts.module.css.js';

// What the Inspector's views share: quiet text tabs, a facts list, colored JSON, the budget meter, a usage line, thin
// bars, and the keys and scrolling of their lists.

export { classes as parts };

export interface TextTab<T extends string> {
  id: T;
  label: string;
  /** A count after the label (Trace 3, Tools 14). */
  count?: MaybeReactive<string | number>;
}

/** Tabs as words: the current one in the text color on the control surface, the others muted. */
export function TextTabs<T extends string>(props: {
  label: string;
  tabs: readonly TextTab<T>[];
  value: Observable<T>;
  onChange(id: T): void;
  class?: string;
}) {
  return (
    <div class={`${classes.tabs} ${props.class ?? ''}`} role="tablist" aria-label={props.label}>
      {props.tabs.map((t) => (
        <button
          type="button"
          role="tab"
          class={classes.tab}
          aria-selected={props.value.pipe(
            map((v) => v === t.id),
            distinctUntilChanged(),
          )}
          onClick={() => props.onChange(t.id)}
        >
          {t.label}
          {t.count !== undefined ? <span class={classes.count}>{t.count}</span> : null}
        </button>
      ))}
    </div>
  );
}

/** A detail's kind (quiet, mono) and where it is (its id), with anything after. */
export function DetailHead(props: { kind: string; where?: string; children?: Child }) {
  return (
    <div class={classes.head}>
      <span class={classes.kind}>{props.kind}</span>
      {props.where ? <span class={classes.where}>{props.where}</span> : null}
      {props.children}
    </div>
  );
}

/** Labelled facts, one per line. */
export function Facts(props: { items: readonly (readonly [string, Child])[] }) {
  return (
    <dl class={classes.facts}>
      {props.items.map(([label, value]) => (
        <div class={classes.fact}>
          <dt>{label}</dt>
          <dd>{value}</dd>
        </div>
      ))}
    </dl>
  );
}

const JsonToken = /("(?:[^"\\]|\\.)*")(\s*:)?|\b(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?|true|false|null)\b/g;

/** JSON, indented, with keys, strings and other values told apart; nothing more. */
export function Json(props: { value: unknown }) {
  const text = JSON.stringify(props.value, null, 2) ?? 'undefined';
  const pre = (<pre class={`${classes.code} ${classes.json}`} />) as HTMLPreElement;
  let at = 0;
  for (const m of text.matchAll(JsonToken)) {
    if (m.index > at) pre.append(text.slice(at, m.index));
    const span = document.createElement('span');
    span.className = m[1] ? (m[2] ? classes.jKey : classes.jStr) : classes.jLit;
    span.textContent = m[1] ?? m[3];
    pre.append(span);
    if (m[2]) pre.append(m[2]);
    at = m.index + m[0].length;
  }
  if (at < text.length) pre.append(text.slice(at));
  return pre;
}

/** A call's or turn's usage in one line (`@molstar/ai/trace` `formatUsage`). */
export function Usage(props: { tokens: TokensView; class?: string }) {
  return <span class={`${classes.usage} ${props.class ?? ''}`}>{formatUsage(props.tokens)}</span>;
}

/** A thin bar: `share` (0..1) of it filled. */
export function Meter(props: { share: number | Observable<number>; class?: string; danger?: boolean }) {
  const width = (s: number) => `${Math.max(0, Math.min(1, s)) * 100}%`;
  const share = props.share;
  return (
    <span class={`${classes.meter} ${props.class ?? ''}`} data-danger={props.danger ? '' : undefined}>
      <span
        class={classes.meterFill}
        style={{ width: typeof share === 'number' ? width(share) : share.pipe(map((s) => width(s))) }}
      />
    </span>
  );
}

/**
 * The latest turn against its budget: the kind of count nearest its ceiling, consumed / ceiling, and its credits.
 * Nothing while there is no trace.
 */
export function BudgetMeter(props: { budget: Observable<BudgetView | undefined> }) {
  return (
    <Show when={props.budget}>
      {(b) => (
        <span
          class={classes.budget}
          data-tip={`${b.scope ? `The ${b.scope}: ` : ''}its budget is counted by kind; this is the kind nearest its ceiling.${b.credits !== undefined ? ' Credits: what it cost, in the app’s credits.' : ''}`}
        >
          {b.ceiling > 0 ? <Meter share={b.consumed / b.ceiling} class={classes.budgetBar} /> : null}
          <span class={classes.budgetText}>
            {b.ceiling > 0 ? `${k(b.consumed)} / ${k(b.ceiling)}` : k(b.consumed)} {b.dimension ?? 'tokens'}
            {b.credits !== undefined ? ` · ${credits(b.credits)}` : ''}
          </span>
        </span>
      )}
    </Show>
  );
}

/** The first line of a text (a result, an entry). */
export const firstLine = (text: string | undefined) => (text ?? '').split('\n', 1)[0] ?? '';

/** The first sentence of a description. */
export function firstSentence(text: string): string {
  const m = /^(.+?[.!?])(\s|$)/.exec(text);
  return (m?.[1] ?? text).replace(/\s+/g, ' ');
}

/**
 * Moves through a list with j/k and the arrows while the page has the keys (not a field, a list or a menu of its own):
 * a listener on the document, so it works whatever in the page has focus.
 */
export function listKeys(root: () => Element | undefined, move: (delta: number) => void) {
  const onKey = (e: KeyboardEvent) => {
    if (e.defaultPrevented || e.altKey || e.ctrlKey || e.metaKey) return;
    const delta = e.key === 'ArrowDown' || e.key === 'j' ? 1 : e.key === 'ArrowUp' || e.key === 'k' ? -1 : 0;
    const el = root();
    if (!delta || !el?.isConnected || !ownsKey(el, e)) return;
    const target = e.target as HTMLElement | null;
    if (
      !el.contains(target) &&
      (target?.isContentEditable ||
        target?.closest?.('input, textarea, select, [role=listbox], [role=menu], [popover]'))
    )
      return;
    // Only from within the page this list is in (or with nothing focused).
    const page = el.closest('[data-msp-page]') ?? el;
    if (target && target !== document.body && !page.contains(target)) return;
    e.preventDefault();
    move(delta);
  };
  on(document, 'keydown', onKey);
}

/**
 * Scrolls a list's viewport (only it: scrollIntoView would move the page too) to show a row, below `top` px of sticky
 * headers; `center` to bring it to the middle.
 */
export function reveal(viewport: HTMLElement | undefined, key: string | undefined, how: 'nearest' | 'center', top = 0) {
  if (!viewport || !key) return;
  let row: HTMLElement | undefined;
  for (const el of viewport.querySelectorAll<HTMLElement>('[data-row]')) {
    if (el.dataset.row === key) {
      row = el;
      break;
    }
  }
  if (!row) return;
  const v = viewport.getBoundingClientRect();
  const r = row.getBoundingClientRect();
  const y = r.top - v.top + viewport.scrollTop;
  if (how === 'center') viewport.scrollTop = y - viewport.clientHeight / 2 + r.height / 2;
  else if (y < viewport.scrollTop + top) viewport.scrollTop = y - top;
  else if (y + r.height > viewport.scrollTop + viewport.clientHeight)
    viewport.scrollTop = y + r.height - viewport.clientHeight;
}
