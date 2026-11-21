/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { extensibilityStory } from './extensibility-root.js';

// Extensibility (Design/Extension openings, part 1): an example extension, `demo`, registered and taken away at run
// time over one spec that names it.

const meta: Meta = { title: 'Features/Extensibility' };
export default meta;

/**
 * The viewer starts without the `demo` extension: the sites track draws its fallback (`ranges`, a bar a site), the
 * domains draw in the layer type's default colors since the `demo:warm` scheme is unknown, and the listing under the
 * canvas says both (`seq.view.layer.extension` once for the extension, `seq.render.theme.scheme`). "Register the
 * extension" adds it at run time: the sites become lollipops (short stems: the section's `demo:stems` key), the domains
 * turn warm, the issues go, and Alt+double-click anywhere runs its command (the sites are shown). "Take the extension
 * away" is its undo: the fallback and the issues come back.
 */
export const ToggledAtRunTime: StoryObj = {
  render: () =>
    extensibilityStory({
      hint: 'A spec that names an extension’s layer type, its key on a section and its scheme, in a viewer without it: the sites are drawn as their fallback (ranges) and the listing says what is missing. Register the extension: the sites become lollipops and the domains warm; Alt+double-click shows the sites. Take it away: the fallback comes back.',
    }),
};
