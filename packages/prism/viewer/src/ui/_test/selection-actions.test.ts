/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { describe, expect, it, vi } from 'vitest';
import type { SelectionAction } from '../../model/selection.js';
import type { ViewerUIContext } from '../../spec.js';
import { extensionUIPort } from '../extensions.js';

// An extension's rows of the selection bar's menu (`ui.selectionActions`): in while it is, each running its command.

describe("an extension's selection actions", () => {
  it('join the menu while the extension is mounted, and run their command', () => {
    const open = vi.fn();
    const actions = observableState<readonly SelectionAction[]>([]);
    const port = extensionUIPort(new ShellModel(), { open } as never, {} as ViewerUIContext, actions);
    const remove = port.register('clipping', {
      selectionActions: [
        {
          id: 'clipping.selection-action.clip',
          title: 'Clip…',
          hint: 'Cut the representations it is in',
          icon: 'crop',
          order: 10,
          command: 'clipping.add',
        },
      ],
    });
    expect(actions.value.map((a) => [a.id, a.label, a.hint, a.order])).toEqual([
      ['clipping.selection-action.clip', 'Clip…', 'Cut the representations it is in', 10],
    ]);
    actions.value[0].run();
    expect(open).toHaveBeenCalledWith('clipping.add');
    remove();
    expect(actions.value).toEqual([]);
  });
});
