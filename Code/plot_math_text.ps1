# Math-style text for the System.Drawing plot scripts.
#
# Dot-source this file, then use Measure-MathText / Draw-MathText in place of
# MeasureString / DrawString for any label that contains notation. Markup:
#   \alpha \beta \theta \Theta \Delta \calB   Greek letters and script B
#   _{...}                                   subscript (smaller, lowered)
#   \bf{...}                                 bold, for matrices such as \bf{B}
# Plain text passes through unchanged. Script files stay ASCII because the
# Greek letters are built from code points here.

# A list of pairs, not a hashtable: PowerShell hashtable keys ignore case, which
# would merge \theta and \Theta. String.Replace below is case-sensitive.
$script:MathSymbols = @(
    @('\alpha', [string][char]0x03B1),
    @('\beta',  [string][char]0x03B2),
    @('\theta', [string][char]0x03B8),
    @('\Theta', [string][char]0x0398),
    @('\Delta', [string][char]0x0394),
    @('\calB',  [string][char]0x212C)
)

$script:MathFormat = [System.Drawing.StringFormat]::GenericTypographic.Clone()
$script:MathFormat.FormatFlags = $script:MathFormat.FormatFlags -bor
    [System.Drawing.StringFormatFlags]::MeasureTrailingSpaces

function ConvertTo-MathSegments([string]$Text) {
    foreach ($pair in $script:MathSymbols) { $Text = $Text.Replace($pair[0], $pair[1]) }
    $segments = @()
    $pattern = '_\{([^}]*)\}|\\bf\{([^}]*)\}'
    $pos = 0
    foreach ($m in [regex]::Matches($Text, $pattern)) {
        if ($m.Index -gt $pos) { $segments += @{ Text = $Text.Substring($pos, $m.Index - $pos); Kind = 'normal' } }
        if ($m.Groups[1].Success) { $segments += @{ Text = $m.Groups[1].Value; Kind = 'sub' } }
        else { $segments += @{ Text = $m.Groups[2].Value; Kind = 'bold' } }
        $pos = $m.Index + $m.Length
    }
    if ($pos -lt $Text.Length) { $segments += @{ Text = $Text.Substring($pos); Kind = 'normal' } }
    return $segments
}

function Get-SegmentFont($Font, [string]$Kind) {
    switch ($Kind) {
        'sub'  { return New-Object System.Drawing.Font $Font.FontFamily, ($Font.Size * 0.68), $Font.Style }
        'bold' { return New-Object System.Drawing.Font $Font.FontFamily, $Font.Size, ($Font.Style -bor [System.Drawing.FontStyle]::Bold) }
        default { return $Font }
    }
}

function Measure-MathText($Graphics, [string]$Text, $Font) {
    $width = 0.0
    foreach ($seg in (ConvertTo-MathSegments $Text)) {
        $f = Get-SegmentFont $Font $seg.Kind
        $width += $Graphics.MeasureString($seg.Text, $f, [System.Drawing.PointF]::Empty, $script:MathFormat).Width
        if ($f -ne $Font) { $f.Dispose() }
    }
    $height = $Graphics.MeasureString("Ag", $Font, [System.Drawing.PointF]::Empty, $script:MathFormat).Height
    return New-Object System.Drawing.SizeF ([float]$width), ([float]$height)
}

function Draw-MathText($Graphics, [string]$Text, $Font, $Brush, [float]$X, [float]$Y) {
    $x = $X
    foreach ($seg in (ConvertTo-MathSegments $Text)) {
        $f = Get-SegmentFont $Font $seg.Kind
        # Subscripts sit lower so their baseline drops below the main baseline.
        $dy = if ($seg.Kind -eq 'sub') { $Font.GetHeight($Graphics) * 0.42 } else { 0 }
        $Graphics.DrawString($seg.Text, $f, $Brush, $x, $Y + $dy, $script:MathFormat)
        $x += $Graphics.MeasureString($seg.Text, $f, [System.Drawing.PointF]::Empty, $script:MathFormat).Width
        if ($f -ne $Font) { $f.Dispose() }
    }
}
