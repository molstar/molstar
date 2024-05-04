# Plugin State Representation

The state of the plugin is represented by a JS Object with these components (described in more detail below):

```ts
interface Snapshot {
    // Snapshot of data state tree
    data?: State.Snapshot,
    // Snapshot of behavior state tree
    behaviour?: State.Snapshot,
    // Snapshot for current animation,
    animation?: PluginAnimationManager.Snapshot,
    // Saved camera positions
    cameraSnapshots?: CameraSnapshotManager.StateSnapshot,
    canvas3d?: {
        // Current camera position
        camera?: Camera.Snapshot,
        // Viewport properties such as background color
        viewport?: Canvas3DProps
    }
}
```

When defining the state object, all components are optional, i.e., it is possible to define just the ``data`` component.

Example state is available [here](./example-state.json). In the plugin, it is possible to create and load these objects using ``Download JSON`` 
and ``Open JSON`` buttons in the ``State Snapshots`` section.

# State Tree

The data and behavior of the plugin is stored in a tree data structure implemented in the ``mol-state`` module. This data structure 
strictly separates the definition of the state with its actual instantiation, similar to the relation of HTML and DOM in web browsers.

The snapshot itself is a JS Object with these components

```ts
interface State.Snapshot {
    tree: StateTree.Serialized
}

interface StateTree.Serialized {
    // Transforms serialized in pre-order
    // The first transform must always be a special "root" node with ref: '-=root=-'
    transforms: [Transform.Serialized, StateObjectCell.State][]
}

interface Transform.Serialized {
    // id of the parent transform
    parent: string,
    // id of the corresponding transformer
    transformer: string,
    // parameters of the transform
    params: any,
    // Properties
    props: Transform.Props,
    // reference to this transform node (a unique string, can be UUID)
    ref: string,
    // version of the node (a unique string, can be UUID)
    version: string
}

interface Transform.Props {
    // tag used in state related operation
    tag?: string
    // is the node visible in the UI
    isGhost?: boolean,
    // is the node bound to its parent? (shown as a single node in the UI)
    isBinding?: boolean
}
```

"Built-in" data state transforms and description of their parameters are defined in ``mol-plugin/state/transforms``. Behavior transforms are defined in ``mol-plugin/behavior``.

# Animation State

Defined by ``CameraSnapshotManager.StateSnapshot`` in ``mol-plugin/state/animation/manager.ts``.

# Canvas3D State

Defined by ``Canvas3DParams`` in ``mol-canvas3d/canvas3d.ts``.

# Camera Snapshots

The camera position (defined in ``mol-canvas3d/camera.ts``) is a plain JS object with the type:

```ts
interface Camera.Snapshot {
    mode: Mode, // = 'perspective' | 'orthographic'

    position: Vec3, // array with [x, y, z]
    // Normalized camera direction
    direction: Vec3, // array with [x, y, z]
    up: Vec3, // array with [x, y, z]
    target: Vec3, // array with [x, y, z]

    near: number,
    far: number,
    fogNear: number,
    fogFar: number,

    fov: number,
    zoom: number
}
```

The ``cameraSnapshots`` component of the state are defined in ``mol-plugin/state/camera.ts``

```js
interface CameraSnapshotManager.StateSnapshot {
    entries: Entry[]
}

interface Entry {
    id: UUID, // or any string
    timestamp: string, // timestamp usually in UTC format
    name?: string, // optional name
    description?: string, // optional description
    snapshot: Camera.Snapshot
}
```