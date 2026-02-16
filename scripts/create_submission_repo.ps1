param(
  [string]$Destination = "..\hcocena-submission"
)

$source = Join-Path $PSScriptRoot "..\hCoCena-r-package"
$target = Resolve-Path -LiteralPath $Destination -ErrorAction SilentlyContinue

if (-not $target) {
  New-Item -ItemType Directory -Path $Destination -Force | Out-Null
  $target = Resolve-Path -LiteralPath $Destination
}

Write-Host "Copying package from $source to $target"
Copy-Item -Path (Join-Path $source "*") -Destination $target -Recurse -Force

Write-Host "Done. Initialize git in the destination if you want a standalone Bioconductor submission repo:"
Write-Host "  cd $target"
Write-Host "  git init"
