Add-Type -AssemblyName System.Drawing
. (Join-Path $PSScriptRoot "..\plot_math_text.ps1")

$root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$tableDir = Join-Path $root "Simulation_Results\FCIR_I\tables"
$plotDir = Join-Path $root "Simulation_Results\FCIR_I\plots"
New-Item -ItemType Directory -Force -Path $plotDir | Out-Null

$rows = Import-Csv (Join-Path $tableDir "Rescale_FCIR_I_adaptive_lasso_selection_wide_sparse_Beta_unpenalized.csv")
$Ns = @(50, 100, 200, 400, 800)
$blue = "#0072B2"
$orange = "#D55E00"

function New-Pen([string]$Color, [bool]$Dashed) {
    $pen = New-Object System.Drawing.Pen ([System.Drawing.ColorTranslator]::FromHtml($Color)), 5
    if ($Dashed) { $pen.DashStyle = [System.Drawing.Drawing2D.DashStyle]::Dash }
    return $pen
}

function Get-MetricValues($Rows, [string]$Matrix, [int]$P, [string]$Metric) {
    return @($Rows |
        Where-Object { $_.Matrix -eq $Matrix -and [int]$_.P -eq $P } |
        Sort-Object { [int]$_.N } |
        ForEach-Object {
            $row = $_
            $nNonzero = [double]$row.n_nonzero
            $nZero = [double]$row.n_total - $nNonzero
            switch ($Metric) {
                "FPR" { ($nZero - [double]$row.'no of correct zeros') / $nZero }
                "TPR" { [double]$row.'no of correct non-zeros' / $nNonzero }
                default { [double]$row.$Metric }
            }
        })
}

function Draw-Panel {
    param(
        $Graphics, [float]$X, [float]$Y, [float]$W, [float]$H,
        [string]$Title, [string]$YLabel, [double]$YMax, $Series
    )
    $axisPen = New-Object System.Drawing.Pen ([System.Drawing.Color]::Black), 2
    $gridPen = New-Object System.Drawing.Pen ([System.Drawing.Color]::LightGray), 1
    $font = New-Object System.Drawing.Font "Arial", 30
    $tickFont = New-Object System.Drawing.Font "Arial", 24
    $titleFont = New-Object System.Drawing.Font "Arial", 34, ([System.Drawing.FontStyle]::Bold)
    $left = $X + 150
    $right = $X + $W - 25
    $top = $Y + 75
    $bottom = $Y + $H - 120

    $titleSize = Measure-MathText $Graphics $Title $titleFont
    Draw-MathText $Graphics $Title $titleFont ([System.Drawing.Brushes]::Black) `
        ($X + ($W - $titleSize.Width) / 2) ($Y + 5)

    for ($i = 0; $i -le 4; $i++) {
        $yy = $bottom - ($bottom - $top) * $i / 4
        $Graphics.DrawLine($gridPen, $left, $yy, $right, $yy)
        $value = $YMax * $i / 4
        $label = if ($YMax -le 1.1) { "{0:0.00}" -f $value } else { "{0:0.0}" -f $value }
        $size = $Graphics.MeasureString($label, $tickFont)
        $Graphics.DrawString($label, $tickFont, [System.Drawing.Brushes]::Black,
            $left - $size.Width - 10, $yy - $size.Height / 2)
    }

    $Graphics.DrawLine($axisPen, $left, $top, $left, $bottom)
    $Graphics.DrawLine($axisPen, $left, $bottom, $right, $bottom)
    for ($i = 0; $i -lt $Ns.Count; $i++) {
        $xx = $left + ($right - $left) * $i / ($Ns.Count - 1)
        $Graphics.DrawLine($axisPen, $xx, $bottom, $xx, $bottom + 7)
        $label = [string]$Ns[$i]
        $size = $Graphics.MeasureString($label, $tickFont)
        $Graphics.DrawString($label, $tickFont, [System.Drawing.Brushes]::Black,
            $xx - $size.Width / 2, $bottom + 10)
    }

    $xLabel = "Number of sites, N"
    $xSize = $Graphics.MeasureString($xLabel, $font)
    $Graphics.DrawString($xLabel, $font, [System.Drawing.Brushes]::Black,
        $X + ($W - $xSize.Width) / 2, $bottom + 65)
    $Graphics.TranslateTransform($X + 22, $Y + $H / 2)
    $Graphics.RotateTransform(-90)
    $YLabel = $YLabel -replace '^F1$', 'F_{1}'
    $ySize = Measure-MathText $Graphics $YLabel $font
    Draw-MathText $Graphics $YLabel $font ([System.Drawing.Brushes]::Black) (-$ySize.Width / 2) 0
    $Graphics.ResetTransform()

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
    $tickFont.Dispose()
    $titleFont.Dispose()
}

function Add-Legend($Graphics, $Items, [float]$Y, [float]$CanvasWidth) {
    $font = New-Object System.Drawing.Font "Arial", 29
    $itemWidth = 330
    $startX = ($CanvasWidth - $itemWidth * $Items.Count) / 2
    for ($i = 0; $i -lt $Items.Count; $i++) {
        $x = $startX + $i * $itemWidth
        $Graphics.DrawLine($Items[$i].Pen, $x, $Y + 16, $x + 65, $Y + 16)
        $Graphics.FillEllipse($Items[$i].Brush, $x + 26, $Y + 9, 14, 14)
        Draw-MathText $Graphics $Items[$i].Label $font ([System.Drawing.Brushes]::Black) ($x + 75) $Y
    }
    $font.Dispose()
}

$seriesSpecs = @(
    @{ Matrix = "Beta_mat"; P = 30; Label = "\calB, P = 30"; Color = $blue; Dashed = $false },
    @{ Matrix = "Beta_mat"; P = 60; Label = "\calB, P = 60"; Color = $blue; Dashed = $true },
    @{ Matrix = "A"; P = 30; Label = "\bf{A}, P = 30"; Color = $orange; Dashed = $false },
    @{ Matrix = "A"; P = 60; Label = "\bf{A}, P = 60"; Color = $orange; Dashed = $true }
)

$bitmap = New-Object System.Drawing.Bitmap 2400, 3000
$graphics = [System.Drawing.Graphics]::FromImage($bitmap)
$graphics.Clear([System.Drawing.Color]::White)
$graphics.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$panels = @(
    @{ X = 30; Y = 100; W = 1170; H = 700; Matrix = "Beta_mat"; Metric = "FPR"; Title = "\calB: FPR"; YMax = 0.50 },
    @{ X = 1200; Y = 100; W = 1170; H = 700; Matrix = "A"; Metric = "FPR"; Title = "A: FPR"; YMax = 0.50 },
    @{ X = 30; Y = 800; W = 1170; H = 700; Matrix = "Beta_mat"; Metric = "TPR"; Title = "\calB: TPR"; YMax = 1.05 },
    @{ X = 1200; Y = 800; W = 1170; H = 700; Matrix = "A"; Metric = "TPR"; Title = "A: TPR"; YMax = 1.05 },
    @{ X = 30; Y = 1500; W = 1170; H = 700; Matrix = "Beta_mat"; Metric = "F1"; Title = "\calB: F_{1} score"; YMax = 1.05 },
    @{ X = 1200; Y = 1500; W = 1170; H = 700; Matrix = "A"; Metric = "F1"; Title = "A: F_{1} score"; YMax = 1.05 },
    @{ X = 30; Y = 2200; W = 1170; H = 700; Matrix = "Beta_mat"; Metric = "RMSE"; Title = "\calB: RMSE"; YMax = 2.30 },
    @{ X = 1200; Y = 2200; W = 1170; H = 700; Matrix = "A"; Metric = "RMSE"; Title = "A: RMSE"; YMax = 0.13 }
)
$legendItems = @()
foreach ($panel in $panels) {
    $series = @()
    foreach ($spec in ($seriesSpecs | Where-Object { $_.Matrix -eq $panel.Matrix })) {
        $pen = New-Pen $spec.Color $spec.Dashed
        $brush = New-Object System.Drawing.SolidBrush ([System.Drawing.ColorTranslator]::FromHtml($spec.Color))
        $item = @{
            Pen = $pen
            Brush = $brush
            Label = $spec.Label
            Values = (Get-MetricValues $rows $spec.Matrix $spec.P $panel.Metric)
        }
        $series += $item
        if ($panels.IndexOf($panel) -lt 2) { $legendItems += $item }
    }
    Draw-Panel $graphics $panel.X $panel.Y $panel.W $panel.H $panel.Title $panel.Metric $panel.YMax $series
}
Add-Legend $graphics $legendItems 20 2400
$bitmap.Save((Join-Path $plotDir "FCIR_I_sparse_selection_metrics.png"), [System.Drawing.Imaging.ImageFormat]::Png)
$graphics.Dispose()
$bitmap.Dispose()

$alphaRows = $rows | Where-Object { $_.Matrix -eq "alpha_0" }
$bitmap = New-Object System.Drawing.Bitmap 1800, 950
$graphics = [System.Drawing.Graphics]::FromImage($bitmap)
$graphics.Clear([System.Drawing.Color]::White)
$graphics.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$series = @()
foreach ($spec in @(
    @{ P = 30; Label = "P = 30"; Color = $blue; Dashed = $false },
    @{ P = 60; Label = "P = 60"; Color = $orange; Dashed = $true }
)) {
    $pen = New-Pen $spec.Color $spec.Dashed
    $brush = New-Object System.Drawing.SolidBrush ([System.Drawing.ColorTranslator]::FromHtml($spec.Color))
    $series += @{
        Pen = $pen
        Brush = $brush
        Label = $spec.Label
        Values = (Get-MetricValues $alphaRows "alpha_0" $spec.P "RMSE")
    }
}
Draw-Panel $graphics 20 80 1760 800 "\alpha_{0}" "RMSE" 0.13 $series
Add-Legend $graphics $series 20 1800
$bitmap.Save((Join-Path $plotDir "FCIR_I_alpha0_RMSE.png"), [System.Drawing.Imaging.ImageFormat]::Png)
$graphics.Dispose()
$bitmap.Dispose()
