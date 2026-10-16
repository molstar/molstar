/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CommandManager, HelpPanel, observableState, render } from '@molstar/prism-ui';
import { afterEach, describe, expect, it } from 'vitest';
import { CanvasBindings } from '../../../ui/bindings.js';
import { HelpProse } from '../help-prose.js';

// Help's prose (parity.md §5): the basics above the bindings.

afterEach(() => {
  document.body.innerHTML = '';
});

describe("Help's prose", () => {
  it('says the basics before the bindings', () => {
    const host = document.createElement('div');
    document.body.append(host);
    // The panel reads the registry's list only.
    const commands = { commands: observableState<readonly never[]>([]).asReadonly } as unknown as CommandManager<any>;
    render(() => <HelpPanel commands={commands} bindings={CanvasBindings} intro={() => <HelpProse />} />, host);
    const titles = [...host.querySelectorAll('dt')].map((d) => d.textContent);
    expect(titles.slice(0, 6)).toEqual(['Load', 'Select', 'Focus', 'Commands', 'Undo', 'Assistant']);
    // The bindings follow.
    expect(titles).toContain('Rotate');
    expect(host.querySelector('[aria-label="Basics"]')?.textContent).toContain('command menu');
  });
});
