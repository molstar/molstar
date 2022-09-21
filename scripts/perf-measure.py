import os
import os.path

found = False

forbidden = [
    'src/cli',
    'src/apps',
    'src/examples',
    'src/servers',
    'src/perf-tests',
    'src/tests',
]

for dirpath, _, filenames in os.walk("../src"):
    dirpath = dirpath.replace('\\', '/')
    if any(f in dirpath for f in forbidden):
        continue

    for filename in [f for f in filenames if f.endswith(".ts") or f.endswith(".tsx")]:
        if '.spec.' in filename or filename.endswith('.d.ts') or filename == "debug.ts":
            continue

        p = os.path.join(dirpath, filename).replace('\\', '/')
        n_slashes = p.count('/') 

        debug_prefix = "../" * (n_slashes - 2)
        debug_import = f"\n\nimport {{ loadCheckpoint }} from '{debug_prefix}mol-util/debug';"
        with open(p, "r", encoding="utf-8") as file:
            text = file.read()

        with open(p, "w", encoding="utf-8") as file:
            text = debug_import + f"\nloadCheckpoint(`{p[7:]}::start`);\n" + text + f"\nloadCheckpoint(`{p[7:]}::end`);\n"
            file.seek(0)
            file.write(text)
            # file.write(debug_import + f"\nloadCheckpoint(`{p[7:]}`);")
