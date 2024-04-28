# Tunnel Visualization Extension
This documentation outlines the usage of the Mol* extension for visualizing tunnels in molecular structures. The extension integrates with Mol* to render 3D representations of tunnels using specified data sources and properties.

## Data Types
The primary data types involved in tunnel visualization are:

### Tunnel
A Tunnel object contains the actual tunnel data necessary for visualization. It consists of:

- `data`: An array of `Profile` objects that describe the tunnel at various points.
- `props`: Properties such as the tunnel's type, ID, and optional labels or descriptions.

### Profile
A `Profile` object in a `Tunnel` holds detailed geometric and physical properties of a tunnel at specific points along its length. These properties include:

- `Charge`: The electric charge at a specific point in the tunnel.
- `Radius`: The overall radius of the tunnel at this point.
- `FreeRadius`: The radius of the tunnel not obstructed by any molecular elements.
- `T`: Temperature factor or a similar property related to the point.
- `Distance`: Distance along the tunnel's path from the start.
- `X`, `Y`, `Z`: Coordinates of the point in 3D space.

These profiles are crucial for understanding the physical and chemical environment inside the tunnel, allowing for detailed analysis and visualization.

## Transformers Usage
The extension uses several transformations to process and visualize tunnel data:

### Tunnels Data Transformer
- `Purpose`: Converts a collection of Tunnel data into a state object.
- `Usage`:
    ```typescript
    update.toRoot().apply(TunnelsDataTransformer, { data: tunnels });
    ```

### Tunnel Data Provider
- `Purpose`: Converts single Tunnel data into a state object for individual processing.
- `Usage`:
    ```typescript
    update.toRoot().apply(TunnelDataTransformer, { 
        data: { 
            data: tunnel.Profile, 
            props: { id: tunnel.Id, type: tunnel.Type } 
        } 
    });
    ```

### Tunnel Shape Provider
- `Purpose`: Provides the shapes for rendering the tunnel based on WebGL context and shape parameters.
- `Usage`:
    ```typescript
    update.apply(TunnelShapeProvider, {
        webgl,
    }).apply(StateTransforms.Representation.ShapeRepresentation3D);
    ```

## Visualization Examples
To help users understand how to use these transformations in practice, include detailed examples:

### Visualizing Multiple Tunnels
This example demonstrates how to visualize multiple tunnels from a fetched dataset.
```typescript
async function runVisualizeTunnels(plugin: PluginContext, url: string = URL) {
    const update = plugin.build();
    const webgl = plugin.canvas3dContext?.webgl;
    const response = await (await fetch(url)).json();
    const tunnels: Tunnel[] = response.Channels.map(channel => ({ data: channel.Profile, props: { id: channel.Id, type: channel.Type } }));

    update.toRoot()
          .apply(TunnelsDataTransformer, { data: tunnels })
          .apply(TunnelsToTunnelTransformer)
          .apply(TunnelShapeProvider, { webgl })
          .apply(StateTransforms.Representation.ShapeRepresentation3D);

    await update.commit();
}
```

### Visualizing a Single Tunnel
This example shows how to visualize a single tunnel.
```typescript
async function runVisualizeTunnel(plugin: PluginContext) {
    const update = plugin.build();
    const webgl = plugin.canvas3dContext?.webgl;
    const response = await (await fetch(URL)).json();
    const tunnel = response.Channels.TransmembranePores_MOLE[0];

    update.toRoot()
          .apply(TunnelDataTransformer, { data: { data: tunnel.Profile, props: { id: tunnel.Id, type: tunnel.Type } } })
          .apply(TunnelShapeProvider, { webgl })
          .apply(StateTransforms.Representation.ShapeRepresentation3D);

    await update.commit();
}
```
