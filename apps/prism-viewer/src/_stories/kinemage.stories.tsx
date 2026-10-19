/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { KinemageExtension } from '@molstar/kinemage-extension/prism-viewer';
import { story } from '@molstar/prism-ui/testing/index';
import { RealViewerStory } from '@molstar/prism-viewer/testing/real-viewer';
import type { PrismViewer } from '@molstar/prism-viewer';
import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';

// The kinemage extension in the Prism viewer. A kinemage (.kin, "kinetic image") is a text file of display lists, not a
// molecule: vectors, balls, dots and triangles in groups and subgroups, masters (named switches over lists across
// groups), views (camera presets) and animate frames (groups shown one at a time). MolProbity, Prekin and Probe write
// them (all-atom contacts, validation outliers); KiNG and Mage view them. The extension gives the plugin the format, so a
// .kin loads like any file (`load.files`, drag and drop, Open), and the `kinemage` API area: its operations are the
// v5 right-panel controls, in Prism the Kinemage panel (a tab on the rail beside Structure while a kinemage is loaded, opened with the file), the
// Kinemage commands of the command menu (⌘K: "kin", "master", "group") and the assistant's tools. The file is data/examples/1crn.kin: crambin (data/examples/1crn.cif) written as Prekin would, with
// groups, subgroups, masters, views, animate frames, vectors, balls and dots. No network.

const meta: Meta = { title: 'Viewer/With extensions/Kinemage', parameters: { layout: 'fullscreen' } };
export default meta;

const KinUrl = new URL('../../../../data/examples/1crn.kin', import.meta.url).href;

/** The viewer's default spec with the kinemage extension: all an app needs to show .kin files. */
const withKinemage = (base: PrismViewerSpec): PrismViewerSpec => ({
  ...base,
  extensions: [...(base.extensions ?? []), KinemageExtension],
});

/** The .kin fetched and opened as a file (what drag and drop does), then the story's steps. */
function kinemageStory(steps?: (v: PrismViewer) => Promise<void>): StoryObj {
  return {
    render: () =>
      story(() => (
        <RealViewerStory
          spec={withKinemage}
          setup={async (v) => {
            const text = await (await fetch(KinUrl)).text();
            const loaded = await v.api.load.files({ files: [new File([text], '1crn.kin')] });
            if (!loaded.ok) throw new Error(loaded.message);
            await steps?.(v);
          }}
        />
      )),
  };
}

const ok = (r: { ok: boolean; message?: string }) => {
  if (!r.ok) throw new Error(r.message);
};

/**
 * Crambin as a kinemage, as it opens: mainchain (white) and sidechain (cyan) vectors, the three disulfides (yellow
 * vectors and balls) with dot shells on their sulfurs, on Kinemage's black background, and the Kinemage panel open on the rail beside Structure.
 * The Structure panel lists it as shapes, not as a structure.
 */
export const Loaded: StoryObj = kinemageStory();

/**
 * The Kinemage panel, drawn as the Structure panel: Next frame in its header, then a tree of its views (a click applies
 * one), its animate frames (the row says which is shown), its masters and its groups (with their subgroups), with eyes. Each control calls a `kinemage` operation as the person: an eye is a step in the
 * history, Undo puts it back. Here the `dots` master is off. Its tab (the kinemage glyph, after Structure) is on the rail
 * only while a kinemage is loaded.
 */
export const Panel: StoryObj = kinemageStory(async (v) => {
  ok(await v.api.view.open({ tab: 'kinemage.tab.kinemage' }));
  ok(await v.api.kinemage!.visibility({ kind: 'master', name: 'dots', visible: false }));
});

/**
 * The Kinemage commands: the command menu open on "Show or hide a kinemage group, subgroup or master" (kind, name,
 * shown). Where the v5 right panel's buttons went; ⌘K and "kin" (or "master", "group") lists the three commands.
 */
export const Commands: StoryObj = kinemageStory(async (v) => {
  v.managers.commands.open('kinemage.visibility');
});

/**
 * Masters: the `sidechain` and `dots` masters off (`kinemage.visibility`), so only the mainchain and the disulfides
 * remain. A master switches the lists that name it, in whichever groups they are.
 */
export const Masters: StoryObj = kinemageStory(async (v) => {
  ok(await v.api.kinemage!.visibility({ kind: 'master', name: 'sidechain', visible: false }));
  ok(await v.api.kinemage!.visibility({ kind: 'master', name: 'dots', visible: false }));
});

/** A view: the kinemage's "Disulfide Cys3-Cys40" (`kinemage.view`), the camera on that bond, orthographic. */
export const View: StoryObj = kinemageStory(async (v) => {
  ok(await v.api.kinemage!.view({ view: 'Disulfide Cys3-Cys40' }));
});

/**
 * Animate: the `1crn` and `CA trace` groups are the kinemage's animate frames, shown one at a time; one step
 * (`kinemage.animate`) shows the next, the CA trace (pink).
 */
export const Animate: StoryObj = kinemageStory(async (v) => {
  ok(await v.api.kinemage!.animate({ sequence: 'animate' }));
});
