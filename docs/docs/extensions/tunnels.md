# Tunnel Visualization Extension
This documentation outlines the usage of the Mol* extension for visualizing tunnels in molecular structures. The extension integrates with Mol* to render 3D representations of tunnels using specified data sources and properties.

The extension is a key component in ChannelsDB (https://channelsdb2.biodata.ceitec.cz/), enabling users to visualize tunnels within molecules directly from the database. While it is used with ChannelsDB, users can also input their own data or connect to different databases, ensuring versatility across various research environments.

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

Example:
```json
"Profile": [
    {
        "Radius": 1.49,
        "FreeRadius": 1.49,
        "T": 0,
        "Distance": 0,
        "X": -19.152,
        "Y": -22.654,
        "Z": -13.034,
        "Charge": 0
    },
    {
        "Radius": 1.524,
        "FreeRadius": 1.524,
        "T": 0.00625,
        "Distance": 0.087,
        "X": -19.162,
        "Y": -22.596,
        "Z": -12.969,
        "Charge": 0
    },
    {
        "Radius": 1.56,
        "FreeRadius": 1.56,
        "T": 0.0125,
        "Distance": 0.174,
        "X": -19.171,
        "Y": -22.539,
        "Z": -12.905,
        "Charge": 0
    }
]
```
## Transformers Usage
The extension uses several transformations to process and visualize tunnel data:

### Tunnels Data Transformer
- `Purpose`: Converts a collection of Tunnel data into a state object.
- `Usage`:
    ```typescript
    update.toRoot().apply(TunnelsFromRawData, { data: tunnels });
    ```

### Tunnel Data Provider
- `Purpose`: Converts single Tunnel data into a state object for individual processing.
- `Usage`:
    ```typescript
    update.toRoot().apply(TunnelFromRawData, {
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
This example ([runVisualizeTunnels](../../../src/extensions/sb-ncbr/tunnels/examples.ts#L19)) demonstrates how to visualize multiple tunnels from a fetched dataset.
```typescript
update.toRoot()
        .apply(TunnelsFromRawData, { data: tunnels })
        .apply(SelectTunnel)
        .apply(TunnelShapeProvider, { webgl })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);
```

### Visualizing a Single Tunnel
This example ([runVisualizeTunnel](../../../src/extensions/sb-ncbr/tunnels/examples.ts#L46)) shows how to visualize a single tunnel.
```typescript
update.toRoot()
        .apply(TunnelFromRawData, {
            data: {
                data: tunnel.Profile,
                props: { id: tunnel.Id, type: tunnel.Type }
            }
        })
        .apply(TunnelShapeProvider, { webgl })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);
```
