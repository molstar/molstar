# wwPDB StructConn extension

The STRUCT_CONN category in the mmCIF file format contains details about the connections between portions of the structure. These can be hydrogen bonds, salt bridges, disulfide bridges and so on (see more at <https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Categories/struct_conn.html>).

**wwPDB StructConn extension** in Mol* provides functionality to retrieve and visualize these connections.

The extension exposes three functions, located in `src/extensions/wwpdb/struct-conn/index.ts`. 

- `getStructConns` - to retrieve struct_conn records from a loaded structure
- `inspectStructConn` - to visualize a struct_conn
- `clearStructConnInspections` - to remove visulizations created by `inspectStructConn`


## Example 1

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
            <button onclick="molstar.PluginExtensions.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'disulf1');">disulf1</button>
            <button onclick="molstar.PluginExtensions.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'disulf2');">disulf2</button>
            <button onclick="molstar.PluginExtensions.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'covale1');">covale1</button>
            <button onclick="molstar.PluginExtensions.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'covale2');">covale2</button>
            <button onclick="molstar.PluginExtensions.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'covale3');">covale3</button>
            <button onclick="molstar.PluginExtensions.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'covale4');">covale4</button>
            <button onclick="molstar.PluginExtensions.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'metalc1');">metalc1</button>
            <button onclick="molstar.PluginExtensions.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'metalc2');">metalc2</button>
            <button onclick="molstar.PluginExtensions.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'metalc3');">metalc3</button>
            <button onclick="molstar.PluginExtensions.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, '5elb', 'metalc4');">metalc4</button>
            <button onclick="molstar.PluginExtensions.wwPDBStructConn.clearStructConnInspections(molstarViewer.plugin, '5elb');">CLEAR</button>
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


## Example 2

This is a more elaborated example, which automatically loads `5elb` (or any PDB entry given in the URL after `?pdb=`), retrieves the list of struct_conns, and creates a button for each struct_conn. 

Be aware that some of the struct_conns may be present in the deposited model but not in the preferred assembly (default view). The presented example will raise a dialog window with error message in such cases, e.g. `disulf6` in entry `5elb`.

```html
<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8" />
        <link rel="icon" href="./favicon.ico" type="image/x-icon">
        <title>Mol* Viewer - StructConn Extension Demo</title>
        <link rel="stylesheet" type="text/css" href="molstar.css" />
    </head>
    <style>
        body { margin: 0px; }
        #app { position: absolute; width: 85%; height: 100%; }
        #controls { position: absolute; right: 0; width: 15%; height: 100%; display: flex; flex-direction: column; overflow-y: scroll; }
        h1 { text-align: center; margin: 12px; font-weight: bold; font-size: 120%; }
        button { margin: 4px; margin-top: 0px; }
    </style>
    <body>
        <div id="app"></div>
        <div id="controls">
            <h1 id="pdb-id">Loading...</h1>
            <button onclick="clearInspections();">CLEAR</button>
        </div>
        <script type="text/javascript" src="./molstar.js"></script>
        <script type="text/javascript">
            var pdbId = window.location.search.match(/[?&]pdb=(\w+)/i)?.[1]?.toLowerCase() ?? '5elb';
            var molstarViewer;
            function inspect(structConnId) {
                if (molstarViewer?.plugin) {
                    molstar.PluginExtensions.wwPDBStructConn.inspectStructConn(molstarViewer.plugin, pdbId, structConnId).then(nSelectedAtoms => {
                        if (nSelectedAtoms < 2) alert('Some of the interacting atoms were not found :(\n(maybe not present in the viewed assembly)');
                    });
                }
            }
            function clearInspections() {
                if (molstarViewer?.plugin) {
                    molstar.PluginExtensions.wwPDBStructConn.clearStructConnInspections(molstarViewer.plugin, pdbId);
                }
            }
            molstar.Viewer.create('app', { layoutIsExpanded: false }).then(viewer => {
                molstarViewer = viewer;
                return viewer.loadPdb(pdbId);
            }).then(() => {
                const structConns = molstar.PluginExtensions.wwPDBStructConn.getStructConns(molstarViewer.plugin, pdbId);
                const controls = document.getElementById('controls');
                for (const structConnId in structConns) {
                    const button = document.createElement('button');
                    button.innerText = structConnId;
                    button.addEventListener('click', () => inspect(structConnId));
                    controls.appendChild(button);
                };
                document.getElementById('pdb-id').innerHTML = pdbId;
            });
        </script>
    </body>
</html>
```
