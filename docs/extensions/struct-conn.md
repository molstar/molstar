# wwPDB StructConn extension

The STRUCT_CONN category in the mmCIF file format contains details about the connections between portions of the structure. These can be hydrogen bonds, salt bridges, disulfide bridges and so on (see more at <https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Categories/struct_conn.html>).

**wwPDB StructConn extension** in Mol* provides functionality to retrieve and visualize these connections.

The extension exposes three functions, located in `src/extensions/wwpdb/struct-conn/index.ts`. 

- `getStructConns` - to retrieve struct_conn records from a loaded structure
- `inspectStructConn` - to visualize a struct_conn
- `clearStructConnInspections` - to remove visulizations created by `inspectStructConn`

The following example is a minimal HTML using this functionality:

```html
<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8" />
        <link rel="icon" href="./favicon.ico" type="image/x-icon">
        <title>Mol* Viewer</title>
        <link rel="stylesheet" type="text/css" href="molstar.css" />
    </head>
    <body style="margin: 0px;">
        <div style="position: absolute; width: 100%; height: 10%; padding-block: 10px;">
            <button onclick="molstar.ExtensionData.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'disulf1');">disulf1</button>
            <button onclick="molstar.ExtensionData.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'disulf2');">disulf2</button>
            <button onclick="molstar.ExtensionData.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'covale1');">covale1</button>
            <button onclick="molstar.ExtensionData.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'covale2');">covale2</button>
            <button onclick="molstar.ExtensionData.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'covale3');">covale3</button>
            <button onclick="molstar.ExtensionData.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'covale4');">covale4</button>
            <button onclick="molstar.ExtensionData.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'metalc1');">metalc1</button>
            <button onclick="molstar.ExtensionData.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'metalc2');">metalc2</button>
            <button onclick="molstar.ExtensionData.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'metalc3');">metalc3</button>
            <button onclick="molstar.ExtensionData.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'metalc4');">metalc4</button>
            <button onclick="molstar.ExtensionData.wwPDBStructConn.clearStructConnInspections(molstarViewer.plugin, '5elb');">CLEAR</button>
        </div>
        <div id="app" style="position: absolute; top: 10%; width: 100%; height: 90%;"></div>
        <script type="text/javascript" src="./molstar.js"></script>
        <script type="text/javascript">
            var molstarViewer;
            molstar.Viewer.create('app', { layoutIsExpanded: false }).then(viewer => {
                molstarViewer = viewer;
                viewer.loadPdb('5elb');
            });
        </script>
    </body>
</html>
```

The PDB ID (`'5elb'`) can be replaced be `undefined`, in which case the functions will apply to the first loaded structure.
