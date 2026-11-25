/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Child, defineElements } from '../dom/index.js';
import { UiElement } from '../elements/ui-element.js';

/** Hosts a story's JSX in a UI root (its styles adopted), so its bindings are disposed when Storybook removes it. */
class StoryRoot extends UiElement<() => Child> {
  static readonly tag = 'msui-story-root';
  protected render(content: () => Child): Child {
    return content();
  }
}

defineElements([StoryRoot]);

/** A Storybook `render` result for JSX content, on the theme's base color. */
export function story(content: () => Child): HTMLElement {
  const el = document.createElement(StoryRoot.tag) as StoryRoot;
  el.style.cssText =
    'display: block; min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
  el.model = content;
  return el;
}

export const page = 'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); max-width: 960px;';
export const heading =
  'margin: var(--msui-space-3) 0 0; font-size: var(--msui-text-md); line-height: var(--msui-leading-md); font-weight: var(--msui-weight-strong);';
export const note = 'margin: 0; color: var(--msui-text-muted);';
export const grid =
  'display: grid; grid-template-columns: repeat(auto-fill, minmax(132px, 1fr)); gap: var(--msui-space-4);';

export function Heading(props: { children: Child }) {
  return <h3 style={heading}>{props.children}</h3>;
}

export function Note(props: { children: Child }) {
  return <p style={note}>{props.children}</p>;
}

/** A color token's swatch, with its name. */
export function Swatch(props: { token: string }) {
  return (
    <div style="display: grid; gap: var(--msui-space-1); font-size: var(--msui-text-xs); line-height: var(--msui-leading-xs);">
      <div
        style={{
          height: '36px',
          borderRadius: 'var(--msui-radius-md)',
          background: `var(${props.token})`,
          boxShadow: 'inset 0 0 0 1px var(--msui-line)',
        }}
      />
      <code style="color: var(--msui-text-muted)">{props.token}</code>
    </div>
  );
}

/** A length token's value as the page resolves it. */
export function resolvedLength(token: string): string {
  const probe = document.createElement('div');
  probe.style.cssText = `position: absolute; visibility: hidden; width: var(${token});`;
  document.body.append(probe);
  const px = `${Math.round(probe.getBoundingClientRect().width * 10) / 10}px`;
  probe.remove();
  return px;
}

/** Rows of tokens with their resolved value, an optional bar of their length, and what they are for. */
export function TokenTable(props: { rows: readonly (readonly [string, string])[]; bars?: boolean }) {
  return (
    <div style="display: grid;">
      {props.rows.map(([token, use]) => (
        <div style="display: grid; grid-template-columns: 200px 56px 136px minmax(0, 1fr); column-gap: var(--msui-space-4); align-items: center; min-height: var(--msui-row);">
          <code style="font-size: var(--msui-text-xs);">{token}</code>
          <span
            class="msui-num"
            ref={(el) => {
              // After it is in the page: the token resolves where it is used.
              queueMicrotask(() => (el.textContent = resolvedLength(token)));
            }}
          />
          <div>
            {props.bars ? (
              <div
                style={{
                  width: `var(${token})`,
                  height: '12px',
                  borderRadius: 'var(--msui-radius-sm)',
                  background: 'var(--msui-accent)',
                }}
              />
            ) : null}
          </div>
          <span style={note}>{use}</span>
        </div>
      ))}
    </div>
  );
}
