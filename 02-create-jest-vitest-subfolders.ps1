# 02-create-jest-vitest-subfolders.ps1
# Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass

Get-ChildItem -Recurse -Directory -Filter '__tests__' | ForEach-Object {
    $jestPath = Join-Path $_.FullName 'jest'
    $vitestPath = Join-Path $_.FullName 'vitest'
    
    New-Item -ItemType Directory -Force -Path $jestPath | Out-Null
    New-Item -ItemType Directory -Force -Path $vitestPath | Out-Null

    Write-Host "Created jest/ and vitest/ inside $($_.FullName)"
}
