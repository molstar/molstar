# Load Trajectory from a Custom Format

This section shows a high level example for loading trajectory from custom data in specialized plugin instances. A more complete solution is available for example in form of the [G3D format extension](https://github.com/molstar/molstar/tree/master/src/extensions/g3d).

## Defining and Using a Custom Transformer

```ts
import { StateTransformer } from 'molstar/lib/mol-state';

const CreateTransformer = StateTransformer.builderFactory('custom-namespace');

export interface CustomTrajectoryData {
    // ...
}

export const TrajectoryFromCustomData = CreateTransformer({
    name: 'trajectory-from-custom-data',
    display: 'Trajectory',
    from: PluginStateObject.Root,
    to: PluginStateObject.Molecule.Trajectory,
    params: {
        data: PD.Value<CustomTrajectoryData>(void 0 as any, { isHidden: true }),
    },
})({
    apply({ params }) {
        return Task.create('Trajectory', async (ctx) => {
            const models = await customParse(params.data, ctx);
            return new PluginStateObject.Molecule.Trajectory(models, {
                label: 'Trajectory',
            });
        });
    },
});
```

The ``customParse`` function can usually be implemented 
by modifying/extending an [existing parser already available in Mol*](https://github.com/molstar/molstar/tree/master/src/mol-model-formats/structure).

To use the transformer:

```ts
const data: CustomTrajectoryData = await (await fetch(url)).json();
const trajectory = await plugin.build().toRoot().apply(TrajectoryFromCustomData, { data }).commit();
// Create the representation
await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
```

## Using Mol* to Download the Data

```ts
export const TrajectoryFromCustomData = CreateTransformer({
    name: 'trajectory-from-custom-data',
    display: 'Trajectory',
    from: PluginStateObject.Data.String, // or PluginStateObject.Data.Binary
    to: PluginStateObject.Molecule.Trajectory,
})({
    apply({ a }) {
        return Task.create('Trajectory', async (ctx) => {
            const models = await customParse(a.data, ctx);
            return new PluginStateObject.Molecule.Trajectory(models, {
                label: 'Trajectory',
            });
        });
    },
});

//////////////

const data = await plugin.builders.data.download({ url, isBinary });
const trajectory = await plugin.build().to(data).apply(TrajectoryFromCustomData, { data }).commit();
await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
```