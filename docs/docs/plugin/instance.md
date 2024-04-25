# Creating Plugin Instance


## Intro

What is a plugin? A plugin is a collection of modules that provide functionality to the `Mol*` UI. The plugin is responsible for managing the state of the viewer, internal and user interactions. It has been a previous point of confusion for new users of `Mol*` to associate the __viewer__ part of the library with what is further referred to as the __plugin__. These two are closely connected in the `molstar-plugin-ui` module, which is the user-facing part of the library and ultimately provides the viewer, but they are ultimately distinct. 


It is recommended that you inspect the general class structure of [`PluginInitWrapper`](https://github.com/molstar/molstar/blob/6edbae80db340134341631f669eec86543a0f1a8/src/mol-plugin-ui/plugin.tsx#L41), [`PluginUIContext`](https://github.com/molstar/molstar/blob/6edbae80db340134341631f669eec86543a0f1a8/src/mol-plugin/context.ts#L71) and [`PluginUIComponent`](https://github.com/molstar/molstar/blob/6edbae80db340134341631f669eec86543a0f1a8/src/mol-plugin-ui/base.tsx#L16) to better understand the flow of data and events in the plugin. 
A passing analogy is that a [ `PluginContext` ](https://github.com/molstar/molstar/blob/6edbae80db340134341631f669eec86543a0f1a8/src/mol-plugin/context.ts#L71) is the engine that powers computation, rendering, events and subscriptions inside the molstar UI. All UI components depend on `PluginContext`. 



There are 4 basic ways of instantiating the Mol* plugin.

## ``Viewer`` wrapper

- The most basic usage is to use the ``Viewer`` wrapper. This is best suited for use cases that do not require much custom behavior and are mostly about just displaying a structure.
- See ``Viewer`` class is defined in [src/apps/viewer/app.ts](https://github.com/molstar/molstar/blob/master/src/apps/viewer/app.ts) for available methods and options.

Example usage without using WebPack:

```HTML
<style>
    #app {
        position: absolute;
        left: 100px;
        top: 100px;
        width: 800px;
        height: 600px;
    }
</style>
<!-- 
    molstar.js and .css are obtained from
    - the folder build/viewer after cloning and building the molstar package 
    - from the build/viewer folder in the Mol* NPM package
-->
<link rel="stylesheet" type="text/css" href="molstar.css" />
<script type="text/javascript" src="./molstar.js"></script>

<div id="app"></div>

<script type="text/javascript">
    molstar.Viewer.create('app', {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: true,
        layoutShowLog: false,
        layoutShowLeftPanel: true,

        viewportShowExpand: true,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,

        pdbProvider: 'rcsb',
        emdbProvider: 'rcsb',
    }).then(viewer => {
        viewer.loadPdb('7bv2');
        viewer.loadEmdb('EMD-30210', { detail: 6 });
    });
</script>
```

When using WebPack (or possibly other build tool) with the Mol* NPM package installed, the viewer class can be imported using 

```ts
import { Viewer } from 'molstar/build/viewer/molstar'

function initViewer(target: string | HTMLElement) {
    return new Viewer(target, { /* options */})
}
```

## ``PluginContext`` with built-in React UI

- For more customization options it is possible to use the [``PluginContext``](https://github.com/molstar/molstar/blob/master/src/mol-plugin/context.ts) directly.
- When creating the plugin instance it is possible to customize the [``PluginSpec``](https://github.com/molstar/molstar/blob/master/src/mol-plugin/spec.ts).
- The default [``PluginSpec``](https://github.com/molstar/molstar/blob/master/src/mol-plugin/spec.ts) is available [here](https://github.com/molstar/molstar/blob/master/src/mol-plugin/spec.ts).
- [``PluginConfig``](https://github.com/molstar/molstar/blob/master/src/mol-plugin/config.ts) object provides additional customization options.
- See the [Viewer State Management](viewer-state.md) section for more information on customizing things like background.
- See the [Data State Management](data-state.md) section for more information on build the state.

```ts
import { DefaultPluginUISpec, PluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import { PluginConfig } from 'molstar/lib/mol-plugin/config';

const MySpec: PluginUISpec = {
    ...DefaultPluginUISpec(),
    config: [
        [PluginConfig.VolumeStreaming.Enabled, false]
    ]
}

async function createPlugin(parent: HTMLElement) {
    const plugin = await createPluginUI({
      target: parent,
      spec: MySpec,
      render: renderReact18
    });

    const data = await plugin.builders.data.download({ url: '...' }, { state: { isGhost: true } });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
    await this.plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');

    return plugin;
}

createPlugin(document.getElementById('app')!); // app is a <div> element with position: relative
```

To use the plugin (with the React UI) inside another React app:

A single-plugin setup is shown the example below. In order to initialize multiple
plugins, each with its own context and viewport, some extra steps are required (docs section to be added).

```ts
import { useEffect, createRef } from "react";
import { createPluginUI } from "molstar/lib/mol-plugin-ui";
import { renderReact18 } from "molstar/lib/mol-plugin-ui/react18";
import { PluginUIContext } from "molstar/lib/mol-plugin-ui/context";
/*  Might require extra configuration,
see https://webpack.js.org/loaders/sass-loader/ for example.
create-react-app should support this natively. */
import "molstar/lib/mol-plugin-ui/skin/light.scss";

declare global {
  interface Window {
    molstar?: PluginUIContext;
  }
}


export function MolStarWrapper() {
  const parent = createRef<HTMLDivElement>();

  // In debug mode of react's strict mode, this code will
  // be called twice in a row, which might result in unexpected behavior.
  useEffect(() => {
    async function init() {
        window.molstar = await createPluginUI({
          target: parent.current as HTMLDivElement,
          render: renderReact18
        });

        const data = await window.molstar.builders.data.download(
          { url: "https://files.rcsb.org/download/3PTB.pdb" }, /* replace with your URL */
          { state: { isGhost: true } }
        );
        const trajectory =
          await window.molstar.builders.structure.parseTrajectory(data, "pdb");
        await window.molstar.builders.structure.hierarchy.applyPreset(
          trajectory,
          "default"
        );
    }
    init();
    return () => {
      window.molstar?.dispose();
      window.molstar = undefined;
    };
  }, []);

  return <div ref={parent} style={{ width: 640, height: 480 }}/>;
}

```


Furthermore, if it is desirable in your project to use the `molstar`'s React UI components, but you wish to alter or rearrange the layout, you should take a look at the signatures of [ `PluginUIComponent` ](https://github.com/molstar/molstar/blob/6edbae80db340134341631f669eec86543a0f1a8/src/mol-plugin-ui/base.tsx#L16) which every "control" subclasses. 


[ `SequenceView` ](https://github.com/molstar/molstar/blob/6edbae80db340134341631f669eec86543a0f1a8/src/mol-plugin-ui/sequence.tsx#L221C4-L221C4), for example, can be used separately from the `PluginUI`. Yet you would need to pass the `PluginUIContext` to it in order for it to observe the changes in the state of the plugin. This can be done via a `PluginContextContainer`:
```typescript
// your_app.plugin: PluginUIContext
...
<div className="your_custom_ui">
  <PluginContextContainer plugin={your_app.plugin}>
    <SequenceView />
  </PluginContextContainer>
</div>
```

## Directly using Mol* React UI

```ts
class MolStarWrapper {
  private resolveInit: () => void;
  initialized = new Promise<boolean>(res => { this.resolveInit = () => res(true); });

  private initCalled = false;
  plugin: PluginUIContext;
  async init() {
    if (this.initCalled) return;
    this.initCalled = true;
    this.plugin = ...;
    await this.plugin.init();
    this.resolveInit();
  }
}

function MolStar({ model }: { model: MolStarWrapper }) {
  const [initialized, setInitialized] = useState(false);
  useEffect(() => {
     async function init() {
       await model.init();
       setInitialized(true);
     }
     init();
  }, [model]);

  if (!initialized) return <>Loading</>;
  return <div style={{ ..., position: 'relative' }}>
    <Plugin plugin={model.plugin} />
  </div>;
}
```

## ``PluginContext`` without built-in React UI

- The [``PluginContext``](https://github.com/molstar/molstar/blob/master/src/mol-plugin/context.ts) can be instantiated without using the default React UI.

```HTML
<div id='molstar-parent' style='position: absolute; top: 0; left: 0; right: 0; bottom: 0'>
    <canvas id='molstar-canvas' style='position: absolute; top: 0; left: 0; right: 0; bottom: 0'></canvas>
</div>
```

```ts
import { DefaultPluginSpec, PluginSpec } from 'molstar/lib/mol-plugin/spec';
import { PluginContext  } from 'molstar/lib/mol-plugin/context';
import { PluginConfig } from 'molstar/lib/mol-plugin/config';

const MySpec: PluginSpec = {
    ...DefaultPluginSpec(),
    config: [
        [PluginConfig.VolumeStreaming.Enabled, false]
    ]
}

async function init() {
    const plugin = new PluginContext(MySpec);
    await plugin.init();

    const canvas = <HTMLCanvasElement> document.getElementById('molstar-canvas');
    const parent = <HTMLDivElement> document.getElementById('molstar-parent');

    if (!plugin.initViewer(canvas, parent)) {
        console.error('Failed to init Mol*');
        return;
    }

    // Example url:"https://files.rcsb.org/download/3j7z.pdb" 
    // Example url:"https://files.rcsb.org/download/5AFI.cif" 
    const data = await plugin.builders.data.download({ url: '...' }, { state: { isGhost: true } });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, format); //format is 'mmcif' or 'pdb' etc.
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
}

```

## ``Canvas3D`` without built-in state management

- The ``PluginContext`` object from the above examples can be completely omitted.
- See [Browser Tests](https://github.com/molstar/molstar/tree/master/src/tests/browser) for example usage.

```ts
const canvas = document.getElementById('canvas'); // parent <canvas> element
const canvas3d = Canvas3D.create(Canvas3DContext.fromCanvas(canvas));
canvas3d.animate();
// use the canvas3d object here
```
