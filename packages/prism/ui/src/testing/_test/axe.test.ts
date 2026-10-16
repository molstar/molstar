// @vitest-environment jsdom

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { axeViolations, expectNoAxeViolations } from '../axe.js';

// The helper itself: that it sees what the views' checks rely on it to see (a check that cannot fail checks nothing).

afterEach(() => {
  document.body.innerHTML = '';
});

function mount(html: string) {
  const el = document.createElement('div');
  el.innerHTML = html;
  document.body.append(el);
  return el;
}

describe('expectNoAxeViolations', () => {
  it('passes a named control and a labelled field', async () => {
    const el = mount('<button type="button" aria-label="Close"></button><label>Name <input /></label>');
    await expectNoAxeViolations(el);
  });

  it('fails a button without a name, an unlabelled field, a bad aria attribute and a duplicate id, with a report', async () => {
    const el = mount(
      '<button type="button"></button><input /><div role="button" aria-pressd="true" tabindex="0">x</div><i id="a"></i><b id="a"></b>',
    );
    const ids = (await axeViolations(el)).map((v) => v.id);
    expect(ids).toContain('button-name');
    expect(ids).toContain('label');
    await expect(expectNoAxeViolations(el)).rejects.toThrow(/button-name/);
  });

  it('sees into a shadow tree (the viewer lives in one)', async () => {
    const host = mount('');
    host.attachShadow({ mode: 'open' }).innerHTML = '<button type="button"></button>';
    expect((await axeViolations(host)).map((v) => v.id)).toContain('button-name');
  });

  it('leaves out the rules jsdom cannot evaluate, and what a view turns off', async () => {
    const el = mount('<button type="button" style="color:#fff;background:#fff">x</button><input />');
    expect((await axeViolations(el)).map((v) => v.id)).not.toContain('color-contrast');
    expect((await axeViolations(el, { disable: ['label'] })).map((v) => v.id)).not.toContain('label');
  });
});
