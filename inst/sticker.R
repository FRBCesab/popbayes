#' 
#' Create an Hexagonal Sticker for the Package
#' 

elephant <- png::readPNG(here::here("inst", "elephant.png"))

p <- ggplot2::ggplot(ggplot2::aes(x = mpg, y = wt), data = mtcars) + 
  ggplot2::geom_point(size = 0.1, colour = "white") + 
  ggplot2::geom_smooth(method = "gam", colour = "white", fill = "white",
                       size = 0.2) + 
  rphylopic::add_phylopic(elephant, alpha = 1, x = 28, y = 4.75, ysize = 6.5, 
                          col = "white") +
  ggplot2::theme_void() + 
  ggpubr::theme_transparent()


hexSticker::sticker(
  
  subplot   = p,
  package   = "popbayes",
  filename  = here::here("man", "figures", "hexsticker.png"),
  dpi       = 600,
  
  p_size    = 6.0,         # Title
  u_size    = 1.2,         # URL
  p_family  = "Aller_Rg",
  
  p_color   = "#FFFFFF",   # Title
  h_fill    = "#2C4370",   # Background
  h_color   = "#5C6F94",   # Border
  u_color   = "#FFFFFF",   # URL
  
  p_x       = 1.00,        # Title
  p_y       = 0.65,        # Title
  s_x       = 1.00,        # Subplot
  s_y       = 1.12,        # Subplot
  
  s_width   = 1.40,        # Subplot
  s_height  = 1.10,        # Subplot
  
  url       = "https://frbcesab.github.io/popbayes",
  
  spotlight = TRUE,
  l_alpha   = 0.10,
  l_width   = 4,
  l_height  = 4
)
