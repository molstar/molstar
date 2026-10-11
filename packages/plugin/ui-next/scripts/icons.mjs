/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Extracts the inner markup of the Lucide icons used by the UI from `lucide-static` into
// `src/components/icon/_generated/lucide.ts`. Only the used icons are copied; the package is a devDependency.
//
//   node scripts/icons.mjs            regenerate
//   node scripts/icons.mjs --check    fail if the generated module is missing or stale

import fs from 'node:fs';
import { createRequire } from 'node:module';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const packageRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const output = path.join(packageRoot, 'src', 'components', 'icon', '_generated', 'lucide.ts');

const HEADER = `/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */`;

const LUCIDE_VERSION = '1.48.0';

/** Lucide icon names (file names in `lucide-static/icons`). */
const ICONS = [
  'house',
  'layers',
  'list-tree',
  'settings',
  'rotate-ccw',
  'camera',
  'save',
  'mouse-pointer-click',
  'maximize',
  'maximize-2',
  'minimize',
  'fullscreen',
  'eye',
  'eye-off',
  'trash-2',
  'chevron-right',
  'chevron-down',
  'chevron-left',
  'chevron-up',
  'chevrons-left',
  'chevrons-right',
  'x',
  'plus',
  'minus',
  'play',
  'clapperboard',
  'square',
  'pause',
  'skip-back',
  'skip-forward',
  'step-back',
  'step-forward',
  'undo-2',
  'redo-2',
  'search',
  'info',
  'triangle-alert',
  'circle-alert',
  'circle-check',
  'loader-circle',
  'ellipsis',
  'ellipsis-vertical',
  'download',
  'upload',
  'folder-open',
  'file',
  'link',
  'clipboard',
  'copy',
  'palette',
  'paintbrush',
  'box',
  'crosshair',
  'focus',
  'monitor-cog',
  'axis-3d',
  'ruler',
  'sun',
  'moon',
  'sparkles',
  'sliders-horizontal',
  'panel-right',
  'panel-left',
  'panel-bottom',
  'panel-top',
  'grip-vertical',
  'pin',
  'history',
  'keyboard',
  'circle-help',
  'image',
  'film',
  'share-2',
  'zoom-in',
  'zoom-out',
  'scan',
  'locate',
  'orbit',
  'rotate-3d',
  'atom',
  'glasses',
  'lightbulb',
  'command',
  'text-cursor-input',
  'list',
  'squares-unite',
  'squares-subtract',
  'squares-intersect',
  'squares-exclude',
  'pencil',
  'mouse-pointer-2',
  'scroll-text',
  'check',
  'hexagon',
  'arrow-up',
  'arrow-down',
  'crop',
  'images',
  'cloud-upload',
  'file-code',
];

const require = createRequire(import.meta.url);
const lucideRoot = path.dirname(require.resolve('lucide-static/package.json'));

/** The child elements of a Lucide SVG, one compact string (`<path d="…"/><circle …/>`). */
function innerMarkup(name, svg) {
  const open = svg.match(/<svg\b[^>]*>/);
  const close = svg.lastIndexOf('</svg>');
  if (!open || close < 0) throw new Error(`${name}: not an SVG`);
  const inner = svg.slice(open.index + open[0].length, close);
  return inner
    .split('\n')
    .map((line) => line.trim())
    .filter(Boolean)
    .join('')
    .replace(/\s+\/>/g, '/>')
    .replace(/\s{2,}/g, ' ');
}

function generate() {
  const version = JSON.parse(fs.readFileSync(path.join(lucideRoot, 'package.json'), 'utf8')).version;
  if (version !== LUCIDE_VERSION) throw new Error(`lucide-static ${LUCIDE_VERSION} expected, found ${version}`);
  const names = [...new Set(ICONS)].sort();
  const lines = [
    HEADER,
    '',
    `// Generated from lucide-static ${LUCIDE_VERSION} (ISC License, Copyright (c) Lucide Contributors) by scripts/icons.mjs; do not edit.`,
    '',
    'export const LucideIcons = {',
  ];
  for (const name of names) {
    const file = path.join(lucideRoot, 'icons', `${name}.svg`);
    if (!fs.existsSync(file)) throw new Error(`Unknown Lucide icon "${name}" (${file})`);
    const markup = innerMarkup(name, fs.readFileSync(file, 'utf8'));
    if (!markup) throw new Error(`${name}: empty icon`);
    lines.push(`  '${name}': '${markup.replace(/\\/g, '\\\\').replace(/'/g, "\\'")}',`);
  }
  lines.push('} as const;', '');
  return lines.join('\n');
}

const content = generate();
if (process.argv.includes('--check')) {
  const current = fs.existsSync(output) ? fs.readFileSync(output, 'utf8') : undefined;
  if (current !== content) {
    console.error(
      `${path.relative(packageRoot, output)} is ${current === undefined ? 'missing' : 'stale'}; run "pnpm icons"`,
    );
    process.exit(1);
  }
} else {
  fs.mkdirSync(path.dirname(output), { recursive: true });
  fs.writeFileSync(output, content);
}
