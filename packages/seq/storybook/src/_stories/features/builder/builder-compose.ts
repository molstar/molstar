/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { printSexpr } from '@molstar/core/repr/codecs/sexpr';
import type { Issue, Node } from '@molstar/core/repr/tree';
import { SeqLanguage } from '@molstar/seqx-core/lang';
import type { ExampleView } from '@molstar/seqx-examples';
import * as asData from './as-data.js';
import asDataSource from './as-data.ts?raw';
import * as check from './check.js';
import checkSource from './check.ts?raw';
import * as firstView from './first-view.js';
import firstViewSource from './first-view.ts?raw';
import * as genome from './genome.js';
import genomeSource from './genome.ts?raw';
import * as handles from './handles.js';
import handlesSource from './handles.ts?raw';
import * as themes from './themes.js';
import themesSource from './themes.ts?raw';

// The Builder stories' cases, DOM-free (the story and its headless twin read them). The code a story shows is the
// module's own source, imported with Vite's `?raw` from the same file the view is built from: one copy, so the code
// shown and the view cannot drift.

/** A panel beside the view: a title and its text (monospace). */
export interface Panel {
  readonly title: string;
  readonly text: string;
}

export interface BuilderCase extends ExampleView {
  /** The module's source, without its license header. */
  readonly source: string;
  /** What is printed beside the view, from the same spec. */
  panels(): Panel[];
}

/** A module's source as shown: its license header left out. */
export const shown = (source: string) => source.replace(/^\/\*\*[\s\S]*?\*\/\s*/, '');

export const formatIssues = (issues: readonly Issue[]): string =>
  issues.length
    ? issues.map((i) => `${i.severity} ${i.code}${i.address ? ` at ${i.address}` : ''}\n  ${i.message}`).join('\n')
    : 'No issues.';

const make = (
  id: string,
  title: string,
  summary: string,
  build: () => Node,
  source: string,
  panels: () => Panel[] = () => [],
): BuilderCase => ({ dataset: 'builder', id, title, summary, build, source: shown(source), panels });

export const BuilderCases: readonly BuilderCase[] = [
  make(
    'first-view',
    'A first view',
    'The smallest spec that draws: a sequence from letters, one track, one layer.',
    firstView.build,
    firstViewSource,
  ),
  make(
    'handles',
    'Handles and references',
    'Every call returns the address of what it adds, typed by its kind; pass it where a reference is wanted (a layer over its data, the data over its sequence). A section handed to a data reference is a type error.',
    handles.build,
    handlesSource,
  ),
  make(
    'themes',
    'Layers and themes',
    'A layer’s channels take themes: a scheme on the letters, a field of the data (categorical) on the ranges, rules over a measure on the values.',
    themes.build,
    themesSource,
  ),
  make(
    'genome',
    'A genome view',
    'An assembly from its sequences, a genome section of it (a block a sequence, a ruler in bases), and one chromosome of it, the locus.',
    genome.build,
    genomeSource,
  ),
  make(
    'check',
    'check()',
    'A mistake the types cannot catch: any string is a reference, so an address written out that names nothing compiles. `check()` says what is wrong, the view shows what is there.',
    check.build,
    checkSource,
    () => [{ title: 'check()', text: formatIssues(check.issues()) }],
  ),
  make(
    'as-data',
    'build() as data',
    'The spec is data: `build()` is a normalized tree (JSON, or the S-expression), the same a host sends without the builder.',
    asData.build,
    asDataSource,
    () => {
      const spec = asData.build();
      return [
        { title: 'JSON', text: JSON.stringify(spec, null, 2) },
        { title: 'S-expression', text: printSexpr(spec, { language: SeqLanguage }) },
      ];
    },
  ),
];
