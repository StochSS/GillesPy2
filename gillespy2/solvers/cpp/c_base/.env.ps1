function New-PathIfNotExists
{
    param(
        [string] $Path
    )

    if ( -Not (Test-Path "${Path}"))
    {
        New-Item -ItemType Directory -Path "$Path"
        return $false
    }
    return $true
}

New-PathIfNotExists -Path "obj"
New-PathIfNotExists -Path "bin"

if ( -Not (Test-Path dev))
{
    Copy-Item -Recurse template -Destination dev
}

$env:OBJ_DIR = "obj"
$env:TEMPLATE_DIR = "dev"
$env:OUTPUT_DIR = "bin"
$env:SOLVER = "hybrid"
