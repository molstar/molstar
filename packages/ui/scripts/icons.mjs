/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Extracts the inner markup of the Lucide icons used by the UI (styles §8) from `lucide-static` into
// `src/components/icon/_generated/lucide/<name>.ts`, a module an icon (`lucide/eye-off.ts` exports `lucideEyeOff`), so a
// consumer that needs a few imports only those modules and a bundler (splitting by module) carries no other; and
// `src/components/icon/_generated/lucide.ts`, the map of them all by name (`LucideIcons`, what `Icon` reads), the same
// markup inline (one module, not 133 imports: what `Icon` loads stays one module). Only the used icons are copied; the
// package is a devDependency.
//
//   node scripts/icons.mjs            regenerate
//   node scripts/icons.mjs --check    fail if the generated module is missing or stale

import fs from 'node:fs';
import { createRequire } from 'node:module';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const packageRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const outDir = path.join(packageRoot, 'src', 'components', 'icon', '_generated');
const output = path.join(outDir, 'lucide.ts');
const iconsDir = path.join(outDir, 'lucide');

const HEADER = fs.readFileSync(fileURLToPath(import.meta.url), 'utf8').split('*/')[0] + '*/';

const LUCIDE_VERSION = '1.48.0';

/** Lucide icon names (file names in `lucide-static/icons`). */
const ICONS = [
  'house',
  'terminal',
  'locate-off',
  'circle',
  'cylinder',
  'cone',
  'flip-vertical-2',
  'scissors',
  'shield',
  'scan',
  'swatch-book',
  'book-open',
  'layers',
  'list-tree',
  'settings',
  'rotate-ccw',
  'camera',
  'save',
  'mouse-pointer-click',
  'maximize',
  'maximize-2',
  'minimize-2',
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
  'volume-2',
  'volume-x',
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
  'message-square',
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
  'sun-moon',
  'ban',
  'align-left',
  'chart-gantt',
  'hourglass',
  'sparkles',
  'shapes',
  'sliders-horizontal',
  'panel-right',
  'panel-left',
  'panel-bottom',
  'panel-top',
  'grip-vertical',
  'pin',
  'cloud',
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
  'glasses',
  'lightbulb',
  'command',
  'text-cursor-input',
  'layout-grid',
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
  'gauge',
  'app-window',
  'toggle-right',
  'wrench',
  'tag',
  'move-3d',
  'refresh-cw',
  'list-restart',
  // History (origins of a step).
  'history',
  'user',
  'bot',
  // Representation types beside the molecular glyphs (a segmented picker): putty, lines.
  'worm',
  'line-squiggle',
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

/** An icon's export name: `lucide` and its name in camel case (`eye-off`: `lucideEyeOff`). */
const exportName = (name) => `lucide${name.replace(/(^|-)(\w)/g, (_, __, c) => c.toUpperCase())}`;

function generate() {
  const version = JSON.parse(fs.readFileSync(path.join(lucideRoot, 'package.json'), 'utf8')).version;
  if (version !== LUCIDE_VERSION) throw new Error(`lucide-static ${LUCIDE_VERSION} expected, found ${version}`);
  const names = [...new Set(ICONS)].sort();
  const generated = `// Generated from lucide-static ${LUCIDE_VERSION} (ISC License, Copyright (c) Lucide Contributors) by scripts/icons.mjs; do not edit.`;
  const outputs = [];
  const entries = [];
  for (const name of names) {
    const file = path.join(lucideRoot, 'icons', `${name}.svg`);
    if (!fs.existsSync(file)) throw new Error(`Unknown Lucide icon "${name}" (${file})`);
    const markup = innerMarkup(name, fs.readFileSync(file, 'utf8'));
    if (!markup) throw new Error(`${name}: empty icon`);
    const value = markup.replace(/\\/g, '\\\\').replace(/'/g, "\\'");
    outputs.push([
      path.join(iconsDir, `${name}.ts`),
      [HEADER, '', generated, '', `export const ${exportName(name)} = '${value}';`, ''].join('\n'),
    ]);
    entries.push(`  '${name}': '${value}',`);
  }
  const map = [HEADER, '', generated, '', 'export const LucideIcons = {', ...entries, '} as const;', ''];
  outputs.push([output, map.join('\n')]);
  return outputs;
}

const outputs = generate();
if (process.argv.includes('--check')) {
  // no module of an icon no longer listed
  const listed = new Set(outputs.map(([file]) => file));
  const stray = fs.existsSync(iconsDir)
    ? fs
        .readdirSync(iconsDir)
        .map((f) => path.join(iconsDir, f))
        .filter((f) => !listed.has(f))
    : [];
  if (stray.length) {
    console.error(`${stray.map((f) => path.relative(packageRoot, f)).join(', ')} not generated; run "pnpm icons"`);
    process.exit(1);
  }
  for (const [file, content] of outputs) {
    const current = fs.existsSync(file) ? fs.readFileSync(file, 'utf8') : undefined;
    if (current !== content) {
      console.error(
        `${path.relative(packageRoot, file)} is ${current === undefined ? 'missing' : 'stale'}; run "pnpm icons"`,
      );
      process.exit(1);
    }
  }
} else {
  fs.rmSync(iconsDir, { recursive: true, force: true });
  fs.mkdirSync(iconsDir, { recursive: true });
  for (const [file, content] of outputs) fs.writeFileSync(file, content);
}
