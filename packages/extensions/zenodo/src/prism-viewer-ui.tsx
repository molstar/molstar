/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, For, Icon, Show, TextInput } from '@molstar/prism-ui';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import type { ZenodoFileInfo, ZenodoRecordInfo } from './prism-viewer.js';
import { ZenodoSection } from './prism-viewer-stub.js';

// The Zenodo section of Home's import slot (parity.md §6): v5's Zenodo Import, flat. A record id, DOI or link; List
// shows the record's files, each with Load where the viewer reads it. It acts only through the person's API
// (`zenodo.list`, `zenodo.load`); their results say themselves in the status column. A trajectory (a model with its
// coordinates) is `zenodo.load`'s `coordinates`, from the command or the agent.

const sizeOf = (n: number) =>
  n >= 1e9 ? `${(n / 1e9).toFixed(1)} GB` : n >= 1e6 ? `${(n / 1e6).toFixed(1)} MB` : `${Math.ceil(n / 1e3)} kB`;

const Row = 'display: flex; align-items: center; gap: var(--msp-space-2); min-width: 0';
// Text that does not fit fades out (reset.css, msp-fade-overflow).
const Name = [
  'flex: 1; min-width: 0; overflow: hidden; white-space: nowrap',
  'mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent)',
  'animation-name: msp-fade-overflow; animation-fill-mode: both; animation-timeline: scroll(self inline)',
].join('; ');
const Quiet = 'color: var(--msp-color-text-muted); font-size: var(--msp-font-size-sm)';

export function zenodoUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: ZenodoSection,
        slot: 'home.import',
        title: 'Zenodo',
        icon: 'cloud',
        keywords: ['zenodo', 'doi', 'record', 'import'],
        render: domView(() => {
          const record = observableState('');
          const busy = observableState(false);
          const listed = observableState<ZenodoRecordInfo | undefined>(undefined);
          const work = async (job: () => Promise<void>) => {
            if (busy.value) return;
            busy.set(true);
            try {
              await job();
            } finally {
              busy.set(false);
            }
          };
          const list = () =>
            work(async () => {
              const r = await ctx.api.zenodo?.list({ record: record.value });
              listed.set(r?.ok ? r.value : undefined);
            });
          const load = (f: ZenodoFileInfo) =>
            work(async () => {
              await ctx.api.zenodo?.load({ record: listed.value?.id ?? record.value, file: f.name });
            });
          return (
            <div style="display: grid; gap: var(--msp-space-2); min-width: 0">
              <div style={Row}>
                <div style="flex: 1; min-width: 0">
                  <TextInput
                    size="sm"
                    value={record.asReadonly}
                    placeholder="Record id, DOI or link"
                    aria-label="Zenodo record"
                    onInput={(v) => record.set(v)}
                    onSubmit={() => void list()}
                  />
                </div>
                <Button size="sm" disabled={busy.asReadonly} onClick={() => void list()}>
                  List files
                </Button>
              </div>
              <Show when={listed.asReadonly}>
                {(r) => (
                  <div style="display: grid; gap: var(--msp-space-1); min-width: 0">
                    <div style={`${Row}; ${Quiet}`}>
                      <span style={Name}>{r.title}</span>
                      <Button
                        size="xs"
                        variant="ghost"
                        icon={<Icon name="x" />}
                        title="Clear"
                        onClick={() => listed.set(undefined)}
                      />
                    </div>
                    <For each={r.files} fallback={() => <span style={Quiet}>No files.</span>}>
                      {(f) => (
                        <div style={Row}>
                          <span style={Name} title={f.name}>
                            {f.name}
                          </span>
                          <span style={Quiet}>{f.kind ? `${f.kind} · ${sizeOf(f.size)}` : sizeOf(f.size)}</span>
                          <Button
                            size="xs"
                            variant="ghost"
                            icon={<Icon name="download" />}
                            title={f.kind ? `Load ${f.name}` : 'The viewer does not read this format'}
                            disabled={!f.kind || f.kind === 'coordinates' || f.kind === 'topology'}
                            onClick={() => void load(f)}
                          />
                        </div>
                      )}
                    </For>
                  </div>
                )}
              </Show>
            </div>
          );
        }),
      },
    ],
  };
}
