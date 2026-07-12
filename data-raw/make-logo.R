# data-raw/make-logo.R
#
# oalasso hex sticker
# 素材: data-raw/oal.png（O が投げ縄になった OAL のワードマーク・紺地）
# 出力: man/figures/logo.png（README・pkgdown 用）
#
# パッケージルートで実行すること: Rscript data-raw/make-logo.R

required_packages <- c("ggplot2", "hexSticker", "magick")

not_installed <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(not_installed) > 0L) {
  install.packages(not_installed)
}

library(hexSticker)
library(ggplot2)
library(magick)

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

background_color <- "#1C2E44"   # 素材画像の紺をサンプリングした値
border_color <- "#EFE9D9"       # 素材の文字色（生成り）に合わせる
name_color <- "#EFE9D9"

# ------------------------------------------------------------------
# 1. hexSticker 自身に「六角形だけ」を描かせるヘルパー
#    （本体レンダリングと形状が正確に一致するマスク・枠線を得る）
# ------------------------------------------------------------------

empty_plot <- ggplot() + theme_void() + theme_transparent()

render_hex_only <- function(fill, colr, h_size) {
  path <- tempfile(fileext = ".png")
  sticker(
    empty_plot, package = "",
    s_x = 1, s_y = 1, s_width = 0.1, s_height = 0.1,
    h_fill = fill, h_color = colr, h_size = h_size,
    white_around_sticker = FALSE,
    filename = path, dpi = 300
  )
  image_read(path)
}

# ------------------------------------------------------------------
# 2. 本体 → 六角形マスクで切り抜き → 枠線
#    地色を素材の紺に合わせているので画像の矩形境界は見えない
# ------------------------------------------------------------------

tmp <- tempfile(fileext = ".png")
sticker(
  "data-raw/oal.png",
  package = "oalasso",
  p_x = 1.00, p_y = 0.55, p_size = 12,
  p_color = name_color, p_family = "Aller_Rg", p_fontface = "bold",
  s_x = 1.00, s_y = 1.10, s_width = 0.85, s_height = 0.85,
  h_fill = background_color, h_color = NA,
  white_around_sticker = FALSE,
  filename = tmp, dpi = 300
)
img <- image_read(tmp)

mask <- render_hex_only("black", NA, 1.2)
mask_gray <- image_negate(image_flatten(image_background(mask, "white")))
mask_gray <- image_convert(mask_gray, type = "grayscale", matte = FALSE)
img <- image_composite(img, mask_gray, operator = "copyopacity")

ring <- render_hex_only(NA, border_color, 1.8)
img <- image_composite(img, ring, operator = "over")

image_write(img, "man/figures/logo.png")
cat("wrote man/figures/logo.png\n")
