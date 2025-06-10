#### CUSTOM THEME ####

theme_bw_alt <- theme(
  text = element_text(size = 20),           # Base text size
  plot.title = element_text(size = 25, face = "bold", hjust = 0),     # Title size
  axis.title.x = element_text(size = 22),   # X-axis title size
  axis.title.y = element_text(size = 22, angle = 90),   # Y-axis title size
  axis.text.x = element_text(size = 18),      # x-axis text size
  axis.text.y = element_text(size = 18),      # y-axis text size
  legend.text = element_text(size = 18),    # Legend text size
  legend.title = element_text(size = 22),   # Legend title size
  strip.background = element_rect(fill = "white", colour = "black"),
  plot.background = element_rect(fill = "white", colour = "white"),
  panel.background = element_rect(fill = "white", colour = NA), 
  panel.border = element_rect(fill = NA, colour = "white"),
  panel.grid.major = element_line(colour = "grey90", size = 0.2),
  panel.grid.minor = element_line(colour = "grey98", size = 0.5))

theme_void_alt <- theme(
  text = element_text(size = 20),           # Base text size
  plot.title = element_text(size = 25, face = "bold", hjust = 0),     # Title size
  axis.title.x = element_blank(),   # X-axis title size
  axis.title.y = element_blank(),   # Y-axis title size
  axis.text = element_blank(),      # Axis text size
  axis.ticks = element_blank(),
  legend.text = element_text(size = 18),    # Legend text size
  legend.title = element_text(size = 22),   # Legend title size
  strip.background = element_rect(fill = "white", colour = "black"),
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white", colour = "white"), 
  panel.border = element_rect(fill = NA, colour = "white"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())
