# Assign custom conformation to a Model

This document shows how to update model conformation dynamically using the `ModelWithCoordinates` transforms. If this does not work well with your particular use case, it is suggested to write a custom version of `ModelWithCoordinates` with similar usage as outlined in this document.

```ts
async function animateFirstXCoordinateExample(plugin: PluginContext, url: string, format: BuiltInTrajectoryFormat) {
    // Load data
    const _data = await plugin.builders.data.download({ url });
    const trajectory = await plugin.builders.structure.parseTrajectory(_data, format);
    const hierarchy = await this.plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
    if (!hierarchy) return;

    // Insert ModelWithCoordinates cell to be updated in the loop bellow
    const coordinatesNode = await plugin.build().to(hierarchy!.model).insert(ModelWithCoordinates).commit();

    const x0 = hierarchy!.model.data!.atomicConformation.x[0];
    let xOffset = 0;
    async function animateFirstXCoord() {
        // Normally, the whole conformation would come from an API/library call, but here we fake it:
        const { x, y, z } = hierarchy!.model.data!.atomicConformation;
        const nextX = [...(x as number[])];
        nextX[0] = x0 + xOffset;
        xOffset += 0.05;
        if (xOffset > 1) xOffset = 0;

        // Construct new coodinate frame from the data and commit the update.
        // Rest of the state tree will reconcile automatically.
        await plugin.build().to(coordinatesNode).update({
            atomicCoordinateFrame: {
                elementCount: x.length,
                time: { value: 0, unit: 'step' },
                xyzOrdering: { isIdentity: true },
                x: nextX,
                y,
                z,
            }
        }).commit();

        requestAnimationFrame(animateFirstXCoord);
    }
    animateFirstXCoord();
}

// animateFirstXCoordinateExample('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/2244/record/SDF/?record_type=3d', 'sdf');
```