# 01-create-__tests__-folders.ps1

Get-ChildItem -Recurse -Filter *.ts -Exclude *.d.ts | ForEach-Object {
    $dir = $_.DirectoryName
    $testDir = Join-Path $dir '__tests__'
    if (-not (Test-Path $testDir)) {
        New-Item -ItemType Directory -Path $testDir | Out-Null
    }
}
