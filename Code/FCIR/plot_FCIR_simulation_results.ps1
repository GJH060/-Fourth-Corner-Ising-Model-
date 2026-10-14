Add-Type -AssemblyName System.Drawing
. (Join-Path $PSScriptRoot "..\plot_math_text.ps1")

$root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$tableDir = Join-Path $root "Simulation_Results\FCIR\tables"
$plotDir = Join-Path $root "Simulation_Results\FCIR\plots"
New-Item -ItemType Directory -Force -Path $plotDir | Out-Null

$selection = Import-Csv (Join-Path $tableDir "adaptive_lasso_selection_wide_unpenalized_Sparse.csv")
$dense = Import-Csv (Join-Path $tableDir "adaptive_lasso_dense_diagnostics_wide_unpenalized_Sparse.csv")
$Ns = @(50, 100, 200, 400, 800)

function New-Pen([string]$color, [float]$width, [bool]$dashed) {
    $pen = New-Object System.Drawing.Pen ([System.Drawing.ColorTranslator]::FromHtml($color)), $width
    if ($dashed) { $pen.DashStyle = [System.Drawing.Drawing2D.DashStyle]::Dash }
    return $pen
}

function Draw-Panel {
    param(
        $Graphics, [float]$X, [float]$Y, [float]$W, [float]$H,
        [string]$Title, [string]$YLabel, [double]$YMax, $Series,
        $ReferenceValue = $null, [string]$ReferenceLabel = ""
    )
    $axisPen = New-Object System.Drawing.Pen ([System.Drawing.Color]::Black), 2
    $gridPen = New-Object System.Drawing.Pen ([System.Drawing.Color]::LightGray), 1
    $font = New-Object System.Drawing.Font "Arial", 32
    $smallFont = New-Object System.Drawing.Font "Arial", 26
    $titleFont = New-Object System.Drawing.Font "Arial", 38, ([System.Drawing.FontStyle]::Bold)
    $left = $X + 155
    $right = $X + $W - 25
    $top = $Y + 85
    $bottom = $Y + $H - 130

    $titleSize = Measure-MathText $Graphics $Title $titleFont
    Draw-MathText $Graphics $Title $titleFont ([System.Drawing.Brushes]::Black) `
        ($X + ($W - $titleSize.Width) / 2) ($Y + 5)

    for ($i = 0; $i -le 4; $i++) {
        $yy = $bottom - ($bottom - $top) * $i / 4
        $Graphics.DrawLine($gridPen, $left, $yy, $right, $yy)
        $value = $YMax * $i / 4
        $label = if ($YMax -le 1.1) { "{0:0.00}" -f $value } else { "{0:0.0}" -f $value }
        $size = $Graphics.MeasureString($label, $smallFont)
        $Graphics.DrawString($label, $smallFont, [System.Drawing.Brushes]::Black,
            $left - $size.Width - 10, $yy - $size.Height / 2)
    }

    $Graphics.DrawLine($axisPen, $left, $top, $left, $bottom)
    $Graphics.DrawLine($axisPen, $left, $bottom, $right, $bottom)

    for ($i = 0; $i -lt $Ns.Count; $i++) {
        $xx = $left + ($right - $left) * $i / ($Ns.Count - 1)
        $Graphics.DrawLine($axisPen, $xx, $bottom, $xx, $bottom + 7)
        $label = [string]$Ns[$i]
        $size = $Graphics.MeasureString($label, $smallFont)
        $Graphics.DrawString($label, $smallFont, [System.Drawing.Brushes]::Black,
            $xx - $size.Width / 2, $bottom + 10)
    }

    $xLabel = "Number of sites, N"
    $xSize = $Graphics.MeasureString($xLabel, $font)
    $Graphics.DrawString($xLabel, $font, [System.Drawing.Brushes]::Black,
        $X + ($W - $xSize.Width) / 2, $bottom + 70)
    $Graphics.TranslateTransform($X + 24, $Y + $H / 2)
    $Graphics.RotateTransform(-90)
    $ySize = Measure-MathText $Graphics $YLabel $font
    Draw-MathText $Graphics $YLabel $font ([System.Drawing.Brushes]::Black) (-$ySize.Width / 2) 0
    $Graphics.ResetTransform()

    if ($null -ne $ReferenceValue) {
        $referencePen = New-Object System.Drawing.Pen ([System.Drawing.Color]::DimGray), 3
        $referencePen.DashStyle = [System.Drawing.Drawing2D.DashStyle]::Dot
        $referenceY = $bottom - ($bottom - $top) * [double]$ReferenceValue / $YMax
        $Graphics.DrawLine($referencePen, $left, $referenceY, $right, $referenceY)
        $labelSize = Measure-MathText $Graphics $ReferenceLabel $smallFont
        $Graphics.FillRectangle(
            [System.Drawing.Brushes]::White,
            $right - $labelSize.Width - 10,
            $referenceY - $labelSize.Height - 2,
            $labelSize.Width + 10,
            $labelSize.Height
        )
        Draw-MathText $Graphics $ReferenceLabel $smallFont ([System.Drawing.Brushes]::DimGray) `
            ($right - $labelSize.Width - 5) ($referenceY - $labelSize.Height - 2)
        $referencePen.Dispose()
    }

    foreach ($series in $Series) {
        $points = @()
        for ($i = 0; $i -lt $Ns.Count; $i++) {
            $xx = $left + ($right - $left) * $i / ($Ns.Count - 1)
            $yy = $bottom - ($bottom - $top) * [double]$series.Values[$i] / $YMax
            $points += New-Object System.Drawing.PointF ([float]$xx), ([float]$yy)
        }
        $Graphics.DrawLines($series.Pen, [System.Drawing.PointF[]]$points)
        foreach ($point in $points) {
            $Graphics.FillEllipse($series.Brush, $point.X - 7, $point.Y - 7, 14, 14)
        }
    }

    $axisPen.Dispose()
    $gridPen.Dispose()
    $font.Dispose()
    $smallFont.Dispose()
    $titleFont.Dispose()
}

function Get-Values($Rows, [string]$Matrix, [int]$P, [string]$Metric) {
    return @($Rows |
        Where-Object { $_.Matrix -eq $Matrix -and [int]$_.P -eq $P } |
        Sort-Object { [int]$_.N } |
        ForEach-Object { [double]$_.$Metric })
}

function Add-Legend($Graphics, $Items, [float]$Y, [float]$CanvasWidth) {
    $font = New-Object System.Drawing.Font "Arial", 30
    $itemWidth = 350
    $startX = ($CanvasWidth - $itemWidth * $Items.Count) / 2
    for ($i = 0; $i -lt $Items.Count; $i++) {
        $x = $startX + $i * $itemWidth
        $Graphics.DrawLine($Items[$i].Pen, $x, $Y + 16, $x + 65, $Y + 16)
        $Graphics.FillEllipse($Items[$i].Brush, $x + 26, $Y + 9, 14, 14)
        Draw-MathText $Graphics $Items[$i].Label $font ([System.Drawing.Brushes]::Black) ($x + 75) $Y
    }
    $font.Dispose()
}

$blue = "#0072B2"
$orange = "#D55E00"
$sparseSeriesSpecs = @(
    @{ Matrix = "B"; P = 30; Label = "\bf{B}, P = 30"; Color = $blue; Dashed = $false },
    @{ Matrix = "B"; P = 60; Label = "\bf{B}, P = 60"; Color = $blue; Dashed = $true },
    @{ Matrix = "A"; P = 30; Label = "\bf{A}, P = 30"; Color = $orange; Dashed = $false },
    @{ Matrix = "A"; P = 60; Label = "\bf{A}, P = 60"; Color = $orange; Dashed = $true }
)

$bitmap = New-Object System.Drawing.Bitmap 2400, 1650
$graphics = [System.Drawing.Graphics]::FromImage($bitmap)
$graphics.Clear([System.Drawing.Color]::White)
$graphics.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$panels = @(
    @{ X = 40; Y = 90; W = 1160; H = 740; Metric = "no of correct zeros"; Title = "Correct zeros"; YMax = 3.2; ReferenceValue = 3; ReferenceLabel = "n_{z} = 3" },
    @{ X = 1200; Y = 90; W = 1160; H = 740; Metric = "no of correct non-zeros"; Title = "Correct nonzeros"; YMax = 3.2; ReferenceValue = 3; ReferenceLabel = "n_{nz} = 3" },
    @{ X = 40; Y = 830; W = 1160; H = 740; Metric = "F1"; Title = "F_{1} score"; YMax = 1.05 },
    @{ X = 1200; Y = 830; W = 1160; H = 740; Metric = "RMSE"; Title = "RMSE"; YMax = 0.25 }
)
$legendItems = @()
foreach ($panel in $panels) {
    $series = @()
    foreach ($spec in $sparseSeriesSpecs) {
        $pen = New-Pen $spec.Color 5 $spec.Dashed
        $brush = New-Object System.Drawing.SolidBrush ([System.Drawing.ColorTranslator]::FromHtml($spec.Color))
        $item = @{ Pen = $pen; Brush = $brush; Label = $spec.Label; Values = (Get-Values $selection $spec.Matrix $spec.P $panel.Metric) }
        $series += $item
        if ($panels.IndexOf($panel) -eq 0) { $legendItems += $item }
    }
    Draw-Panel $graphics $panel.X $panel.Y $panel.W $panel.H $panel.Title $panel.Title $panel.YMax $series $panel.ReferenceValue $panel.ReferenceLabel
}
Add-Legend $graphics $legendItems 20 2400
$bitmap.Save((Join-Path $plotDir "FCIR_sparse_selection_metrics.png"), [System.Drawing.Imaging.ImageFormat]::Png)
$graphics.Dispose()
$bitmap.Dispose()

$bitmap = New-Object System.Drawing.Bitmap 2200, 900
$graphics = [System.Drawing.Graphics]::FromImage($bitmap)
$graphics.Clear([System.Drawing.Color]::White)
$graphics.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$denseSpecs = @(
    @{ P = 30; Label = "P = 30"; Color = $blue; Dashed = $false },
    @{ P = 60; Label = "P = 60"; Color = $orange; Dashed = $true }
)
$densePanels = @(
    @{ X = 40; Y = 90; W = 1060; H = 740; Matrix = "beta_0"; Title = "\beta_{0}"; YMax = 0.65 },
    @{ X = 1100; Y = 90; W = 1060; H = 740; Matrix = "alpha_0"; Title = "\alpha_{0}"; YMax = 0.12 }
)
$legendItems = @()
foreach ($panel in $densePanels) {
    $series = @()
    foreach ($spec in $denseSpecs) {
        $pen = New-Pen $spec.Color 5 $spec.Dashed
        $brush = New-Object System.Drawing.SolidBrush ([System.Drawing.ColorTranslator]::FromHtml($spec.Color))
        $item = @{ Pen = $pen; Brush = $brush; Label = $spec.Label; Values = (Get-Values $dense $panel.Matrix $spec.P "RMSE") }
        $series += $item
        if ($densePanels.IndexOf($panel) -eq 0) { $legendItems += $item }
    }
    Draw-Panel $graphics $panel.X $panel.Y $panel.W $panel.H $panel.Title "RMSE" $panel.YMax $series
}
Add-Legend $graphics $legendItems 20 2200
$bitmap.Save((Join-Path $plotDir "FCIR_dense_RMSE.png"), [System.Drawing.Imaging.ImageFormat]::Png)
$graphics.Dispose()
$bitmap.Dispose()
