/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState } from '@molstar/prism-ui/core';
import { render } from '@molstar/prism-ui/dom/render';
import { afterEach, describe, expect, it } from 'vitest';
import type { AIProvider, AISettingsModel, AISetup as Setup } from '../../model.js';
import { AISetup } from '../ai-setup.js';

type Added = Parameters<AISettingsModel['addProvider']>[0];

function mount() {
  const added: Added[] = [];
  const settings: AISettingsModel = {
    providers: observableState<readonly AIProvider[]>([]).asReadonly,
    setup: observableState<Setup | undefined>(undefined).asReadonly,
    setKey: () => undefined,
    removeKey: () => undefined,
    choose: () => undefined,
    refreshModels: async () => undefined,
    listModels: async () => [],
    addProvider: (p) => void added.push(p),
  };
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <AISetup settings={settings} />, container);
  const button = (text: string) =>
    [...container.querySelectorAll('button')].find((b) => b.textContent?.trim() === text) as HTMLButtonElement;
  button('Your own endpoint').click();
  const input = (label: string) => container.querySelector<HTMLInputElement>(`input[aria-label="${label}"]`)!;
  const type = (label: string, value: string) => {
    const el = input(label);
    el.value = value;
    el.dispatchEvent(new Event('input', { bubbles: true }));
  };
  type('Name', 'Lab');
  type('Endpoint', 'http://localhost:11434/v1');
  type('Models', 'llama');
  return { container, added, input, type, add: () => button('Add').click() };
}

afterEach(() => {
  document.body.innerHTML = '';
});

describe('your own endpoint’s limits', () => {
  it('are optional: the fields say the defaults, and an endpoint added without them has none of its own', () => {
    const { input, added, add } = mount();
    expect(input('Context window').placeholder).toContain('32,768');
    expect(input('Max output').placeholder).toContain('4,096');
    add();
    expect(added).toHaveLength(1);
    expect(added[0]).not.toHaveProperty('contextWindow');
    expect(added[0]).not.toHaveProperty('maxOutputTokens');
  });

  it('are kept with the endpoint when given', () => {
    const { type, added, add } = mount();
    type('Context window', '131,072');
    type('Max output', '8192');
    add();
    expect(added[0]).toMatchObject({ label: 'Lab', models: ['llama'], contextWindow: 131_072, maxOutputTokens: 8192 });
  });

  it('are refused unless whole numbers above zero, with the reason shown', () => {
    const { type, added, add, container } = mount();
    type('Context window', 'lots');
    add();
    expect(added).toHaveLength(0);
    expect(container.textContent).toContain('whole numbers of tokens');
  });
});
