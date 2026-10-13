/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The AI packages' tests run apart from the rest (`pnpm jest:ai`): the AI SDK is ESM only, so its JavaScript is
// transformed to CommonJS here, which the rest of the repository does not need. Goes away when Jest moves to ESM
// (.v6/plans/checklist.md).
const base = require('../../package.json').jest;

const esmOnly = ['ai', '@ai-sdk', '@workflow', '@opentelemetry', 'eventsource-parser', '@standard-schema', '@vercel'];

module.exports = {
  ...base,
  rootDir: '../..',
  roots: ['<rootDir>/packages/ai'],
  testPathIgnorePatterns: ['/node_modules/'],
  transform: {
    ...base.transform,
    '\\.m?js$': ['esbuild-jest-transform', { format: 'cjs' }],
  },
  transformIgnorePatterns: [`/node_modules/(?!(\\.pnpm/[^/]+/node_modules/)?(${esmOnly.join('|')})(/|$))`],
};
