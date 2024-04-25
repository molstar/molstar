# Plugin State Server

Provides a simple backend for online storing and sharing of Mol* sessions used by 
[``mol-plugin``](https://github.com/molstar/molstar/tree/master/src/mol-plugin) and 
[``mol-state``](https://github.com/molstar/molstar/tree/master/src/mol-state) modules.

## Example
```sh
node lib/commonjs/servers/plugin-state --workding-folder ~
```

## Usage
| Argument | Description |
| --- | --- |
| `--working-folder` | Working folder path |
| `--port` | Server port. Alternatively, use ENV variable PORT. |
| `--api-prefix` | Server API prefix |
| `--max-states` | Maximum number of states to save |

```sh
node lib/commonjs/servers/plugin-state [-h] --working-folder WORKING_FOLDER [--port PORT] [--api-prefix API_PREFIX] [--max-states MAX_STATES]
```