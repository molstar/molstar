/** @jest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { LucideIcons } from '../_generated/lucide.js';
import { Icon, type IconName } from '../icon.js';
import { MolecularIcons } from '../molecular.js';

const SVG_NS = 'http://www.w3.org/2000/svg';

describe('Icon', () => {
  it('renders an SVG in the SVG namespace with children', () => {
    const svg = Icon({ name: 'house' });
    expect(svg.namespaceURI).toBe(SVG_NS);
    expect(svg.localName).toBe('svg');
    expect(svg.getAttribute('viewBox')).toBe('0 0 24 24');
    expect(svg.getAttribute('stroke')).toBe('currentColor');
    expect(svg.getAttribute('fill')).toBe('none');
    expect(svg.children.length).toBeGreaterThan(0);
    for (const child of Array.from(svg.children)) expect(child.namespaceURI).toBe(SVG_NS);
  });

  it('renders every icon', () => {
    for (const name of [...Object.keys(LucideIcons), ...Object.keys(MolecularIcons)] as IconName[]) {
      expect(Icon({ name }).children.length).toBeGreaterThan(0);
    }
  });

  it('has no name clashes between Lucide and molecular icons', () => {
    const lucide = new Set(Object.keys(LucideIcons));
    for (const name of Object.keys(MolecularIcons)) expect(lucide.has(name)).toBe(false);
  });

  it('applies size and class', () => {
    const def = Icon({ name: 'x' });
    expect(def.getAttribute('width')).toBe('18');
    expect(def.getAttribute('height')).toBe('18');
    const big = Icon({ name: 'x', size: 32, class: 'custom' });
    expect(big.getAttribute('width')).toBe('32');
    expect(big.getAttribute('height')).toBe('32');
    expect(big.getAttribute('class')).toContain('custom');
  });

  it('returns a new node per call that shares structure', () => {
    const a = Icon({ name: 'layers', size: 16 });
    const b = Icon({ name: 'layers', size: 24 });
    expect(a).not.toBe(b);
    expect(a.firstElementChild).not.toBe(b.firstElementChild);
    expect(a.innerHTML).toBe(b.innerHTML);
    expect(a.getAttribute('width')).toBe('16');
    expect(b.getAttribute('width')).toBe('24');
  });

  it('does not leak per-instance attributes into the template', () => {
    Icon({ name: 'search', label: 'Find', class: 'one', size: 40 });
    const next = Icon({ name: 'search' });
    expect(next.getAttribute('aria-label')).toBeNull();
    expect(next.getAttribute('role')).toBeNull();
    expect(next.getAttribute('width')).toBe('18');
    expect(next.getAttribute('class')).not.toContain('one');
  });

  it('is decorative without a label', () => {
    const svg = Icon({ name: 'eye' });
    expect(svg.getAttribute('aria-hidden')).toBe('true');
    expect(svg.hasAttribute('role')).toBe(false);
    expect(svg.hasAttribute('aria-label')).toBe(false);
  });

  it('is an image with a label', () => {
    const svg = Icon({ name: 'eye', label: 'Visible' });
    expect(svg.getAttribute('role')).toBe('img');
    expect(svg.getAttribute('aria-label')).toBe('Visible');
    expect(svg.hasAttribute('aria-hidden')).toBe(false);
  });

  it('throws for unknown names', () => {
    expect(() => Icon({ name: 'nope' as IconName })).toThrow();
  });
});
